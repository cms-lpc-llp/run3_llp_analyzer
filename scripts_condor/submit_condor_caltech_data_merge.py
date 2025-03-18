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
executable = "runAnalyzer.sh"
analyzer = 'llp_MuonSystem_CA_merge'
filesPerJob = 10
ntupler_version = 'V1p19/'
merged_version = "v2"

analyzer_version = 'v30'
HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"


for year in ['2022','2023','2024']:
    if year == "2024":merged_version = "v3"
    else: merged_version = "v2"
    datasetListDir = Analyzer_DIR + f"lists/MergedNtuples/Run3/{ntupler_version}/{year}/{merged_version}/"
    outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/{0}/Data{1}/{2}/".format(ntupler_version, year,analyzer_version)
    #### Run on all signal samples ####
    ## year/isData
    datasetList = OrderedDict()
    samples = os.listdir(datasetListDir)
    for s in samples:
        if "MET" in s:continue
        if not "2024G" in s:continue
        datasetList[s.replace('.txt', '')] = [year, "yes"]
    ############
    
    
    
    for sample in datasetList.keys():
    
        print("Preparing analyzer workflow for dataset: " + sample + "\n")
    
        inputfilelist  = datasetListDir + sample +'.txt'
        if not os.path.exists(inputfilelist):
            print("listfile: " + inputfilelist + " does not exist. skipping.")
            continue
        print(inputfilelist)
        cmd = ["awk", "END{print NR}",inputfilelist]
        nfiles = int(subprocess.check_output(cmd).decode("utf-8"))
        maxjob=int(nfiles/filesPerJob)+1
        mod=int(nfiles%filesPerJob)
        if mod == 0: maxjob -= 1
        outputDirectory = outputDirectoryBase + sample + "/"
        if year == "2024": analyzerTag = "Summer24" #doesn't really matter for data
        elif year == '2023':
            if '2023D' in sample: analyzerTag = "Summer23BPix"
            elif '2023B' in sample or '2023C' in sample: analyzerTag = "Summer23"
        elif year == '2022':
            if '2022B' in sample or '2022C' in sample or '2022D' in sample: analyzerTag = "Summer22"
            elif '2022E' in sample or '2022F' in sample or '2022G' in sample: analyzerTag = "Summer22EE"
        else:
            sys.exit("ANALYSIS TAG NOT FOUND")
        print(year, analyzerTag, sample)
        year = datasetList[sample][0]
        isData = datasetList[sample][1]
        print(maxjob)
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
    
        tmpCondorJDLFile.write("+JobQueue=\"Short\" \n")
        tmpCondorJDLFile.write("RequestMemory = 2000 \n")
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
