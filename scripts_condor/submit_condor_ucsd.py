#!/usr/bin/env python3

import os
import subprocess
from collections import OrderedDict


os.system("mkdir -p submit")
os.system("mkdir -p log")
executable = "runAnalyzer_ucsd.sh"
# TODO(uaf): adjust runAnalyzer_ucsd.sh if UAF-specific changes are needed (SCRAM_ARCH, xrdcp, etc).
analyzer = 'llp_MuonSystem_CA_mdsnano'
filesPerJob = 1
ntupler_version = 'V1p19/Data2023/'
#ntupler_version = "V1p19/MC_Summer22EE/v1/sixie/"


outputDirectoryBase="/ceph/cms/store/group/mds-ml/run3_analyzer_output/"
# TODO(uaf): if you want ntupler/analyzer subfolders, append them here.
HOME = os.getenv('HOME')
X509_USER_PROXY = os.getenv("X509_USER_PROXY", "/tmp/x509up_u81902")
# TODO(uaf): make proxy path fully agnostic if needed.
CMSSW_BASE = "/ceph/cms/store/group/mds-ml/cmssw/CMSSW_14_1_0_pre4"
Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"
# datasetListDir = Analyzer_DIR + "lists/MDSNano/v2/MC_Summer24/"
datasetListFile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "MC_Summer24_test.txt")



#### Run on all signal samples ####
## isData
datasetList = OrderedDict()
sample = os.path.basename(datasetListFile).replace(".txt", "")
# TODO(uaf): switch back to directory scanning when done testing.
datasetList[sample] = "no"
############



for sample in datasetList.keys():

    print("Preparing analyzer workflow for dataset: " + sample + "\n")

    inputfilelist  = datasetListFile
    if not os.path.exists(inputfilelist):
        print("listfile: " + inputfilelist + " does not exist. skipping.")
        continue
    cmd = ["awk", "END{print NR}",inputfilelist]
    nfiles = int(subprocess.check_output(cmd).decode("utf-8"))
    maxjob=int(nfiles/filesPerJob)+1
    mod=int(nfiles%filesPerJob)
    if mod == 0: maxjob -= 1
    # outputDirectory = "/home/users/aaportel/run3_analyzer_output"
    # outputDirectory = "/ceph/cms/store/user/aaportel/run3_analyzer_output"
    outputDirectory = outputDirectoryBase
    analyzerTag = "Summer24"
    isData = datasetList[sample]
    print(f"inputfilelist={inputfilelist}")
    print(f"nfiles={nfiles} filesPerJob={filesPerJob} maxjob={maxjob}")
    print(f"outputDirectory={outputDirectory}")
    print(f"CMSSW_BASE={CMSSW_BASE}")

    #####################################
    #Create Condor JDL file
    #####################################
    os.system("rm -f submit/{}_{}_Job*.jdl".format(analyzer, sample))
    os.system("rm -f log/{}_{}_Job*".format(analyzer, sample))


    jdl_file="submit/{}_{}_{}.jdl".format(analyzer, sample, maxjob)

    tmpCondorJDLFile = open(jdl_file,"w")

    tmpCondorJDLFile.write("Universe = vanilla \n")
    tmpCondorJDLFile.write("Executable = {} \n".format(executable))
    tmpCondorJDLFile.write("Arguments = {} {} {} {} $(ProcId) {} {} {} {} {}/ \n"\
                            .format(analyzer, inputfilelist, isData, filesPerJob, maxjob, outputDirectory, analyzerTag, CMSSW_BASE, HOME))

    tmpCondorJDLFile.write("Log = log/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).log \n".format(analyzer, sample, maxjob))
    tmpCondorJDLFile.write("Output = log/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).out \n".format(analyzer, sample, maxjob))
    tmpCondorJDLFile.write("Error = log/{}_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).err \n".format(analyzer, sample, maxjob))

    # TODO(uaf): confirm desired JobFlavour (espresso/longlunch/workday/etc). Using longlunch for now.
    tmpCondorJDLFile.write("+JobFlavour = \"espresso\" \n")
    tmpCondorJDLFile.write("+DESIRED_Sites = \"T2_US_UCSD\" \n")
    tmpCondorJDLFile.write("RequestMemory = 1000 \n")
    tmpCondorJDLFile.write("RequestCpus = 1 \n")
    # HTCondor units are KB. 2 GB is a safer minimum for CMSSW jobs.
    tmpCondorJDLFile.write("RequestDisk = 2000000 \n")

    tmpCondorJDLFile.write("+RunAsOwner = True \n")
    # tmpCondorJDLFile.write("+InteractiveUser = true \n")
    # TODO(uaf): confirm OS/container choice (rhel7 vs rhel9/alma8) for UAF.
    tmpCondorJDLFile.write("+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel8\" \n")
    tmpCondorJDLFile.write('+WantOS = "el8" \n')
    tmpCondorJDLFile.write('+SingularityBindCVMFS = True \n')
    # tmpCondorJDLFile.write("run_as_owner = True \n")
    tmpCondorJDLFile.write("use_x509userproxy = True \n")
    tmpCondorJDLFile.write("x509userproxy = {} \n".format(X509_USER_PROXY))
    tmpCondorJDLFile.write("should_transfer_files = YES \n")
    tmpCondorJDLFile.write("when_to_transfer_output = ON_EXIT \n")
    tmpCondorJDLFile.write("StreamOut = True \n")
    tmpCondorJDLFile.write("StreamErr = True \n")
    tmpCondorJDLFile.write('Environment = "RUN_DIAG=1" \n')

    tmpCondorJDLFile.write("Queue {} \n".format(maxjob))
    #tmpCondorJDLFile.write("Queue {} \n".format(2))
    tmpCondorJDLFile.close()

    os.system("condor_submit {} --batch-name {}".format(jdl_file, sample))
