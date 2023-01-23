#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict


os.system("mkdir -p submit")
os.system("mkdir -p log")
executable = "run_printJson.sh"
filesPerJob = 1
ntupler_version = 'V1p19/Data2022/'
analyzer_version = 'v1'
outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/displacedJetMuonNtuple_json/Run3/{0}/{1}/".format(ntupler_version, analyzer_version)
HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"
datasetListDir = Analyzer_DIR + "lists/displacedJetMuonNtuple/{}/".format(ntupler_version)


#### Run on all signal samples ####
datasetList = OrderedDict()
samples = os.listdir(datasetListDir)
for sample_temp in samples:
    sample = sample_temp.replace('.txt', '')
    print("Preparing analyzer workflow for dataset: " + sample + "\n")

    inputfilelist  = datasetListDir + sample + '.txt'
    if not os.path.exists(inputfilelist):
        print("listfile: " + inputfilelist + " does not exist. skipping.")
        continue
    cmd = ["awk", "END{print NR}",inputfilelist]
    nfiles = int(subprocess.check_output(cmd).decode("utf-8"))
    maxjob=int(nfiles/filesPerJob)+1
    mod=int(nfiles%filesPerJob)
    if mod == 0: maxjob -= 1
    outputDirectory = outputDirectoryBase + sample + "/"

    #####################################
    #Create Condor JDL file
    #####################################
    os.system("rm -f submit/printJson_{}_Job*.jdl".format(sample))
    os.system("rm -f log/printJson_{}_Job*".format(sample))


    jdl_file="submit/printJson_{}_{}.jdl".format(sample, maxjob)

    tmpCondorJDLFile = open(jdl_file,"w")

    tmpCondorJDLFile.write("Universe = vanilla \n")
    tmpCondorJDLFile.write("Executable = {} \n".format(executable))
    tmpCondorJDLFile.write("Arguments = {} {} $(ProcId) {} {} {} {}/ \n"\
                            .format(inputfilelist, filesPerJob, maxjob, outputDirectory,  CMSSW_BASE, HOME))

    tmpCondorJDLFile.write("Log = log/printJson_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).log \n".format(sample, maxjob))
    tmpCondorJDLFile.write("Output = log/printJson_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).out \n".format(sample, maxjob))
    tmpCondorJDLFile.write("Error = log/printJson_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).err \n".format(sample, maxjob))

    tmpCondorJDLFile.write("+JobQueue=\"Short\" \n")
    tmpCondorJDLFile.write("RequestMemory = 4000 \n")
    tmpCondorJDLFile.write("RequestCpus = 1 \n")
    tmpCondorJDLFile.write("RequestDisk = 4 \n")

    tmpCondorJDLFile.write("+RunAsOwner = True \n")
    tmpCondorJDLFile.write("+InteractiveUser = true \n")
    tmpCondorJDLFile.write("+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7\" \n")
    tmpCondorJDLFile.write('+SingularityBindCVMFS = True \n')
    tmpCondorJDLFile.write("run_as_owner = True \n")
    tmpCondorJDLFile.write("x509userproxy = {}/x509_proxy \n".format(HOME))
    tmpCondorJDLFile.write("should_transfer_files = YES \n")
    tmpCondorJDLFile.write("when_to_transfer_output = ON_EXIT \n")
    tmpCondorJDLFile.write("Queue {} \n".format(maxjob))
    tmpCondorJDLFile.close()

    os.system("condor_submit {} --batch-name {}".format(jdl_file, sample))
