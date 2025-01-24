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
executable = "runMergeNtuple_cache.sh"
filesPerJob = 1

list_path="/storage/af/user/christiw/login-1/christiw/LLP/Run3/CMSSW_10_6_30/src/run3_llp_analyzer/scripts_local/merge_ntuples_commands/"

outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/MergedNtuples/Run3/V1p19/"
HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"


datasetList = OrderedDict() #list of ntuples and point to the correct nano list file
samples = os.listdir(list_path)
print(samples)
for sample in samples:
    if "lists" in sample:continue # run over JETMET dataset
    if not "Muon0-Run2024H-EXOCSCCluster-PromptReco-v1_v1_v2" in sample and not "Muon1-Run2024I-EXOCSCCluster-PromptReco-v2_v1_v1" in sample:continue
    #if not "Muon0-EXOCSCCluster_Run2023C-PromptReco-v1" in sample:continue
    # if not "ggH" in sample:continue
    print("Preparing workflow for dataset: " + sample + "\n")

    inputfilelist  = list_path + sample 

    # check if list files exist
    if not os.path.exists(inputfilelist):
        print("ntuple list file: " + inputfilelist + " does not exist. skipping.")
        continue
    cmd = ["awk", "END{print NR}",inputfilelist]
    nfiles = int(subprocess.check_output(cmd).decode("utf-8"))
    maxjob=int(nfiles/filesPerJob)+1
    mod=int(nfiles%filesPerJob)
    if mod == 0: maxjob -= 1

    year = sample[sample.find("Run")+3:sample.find("Run")+7]
    outputDirectory = f"{outputDirectoryBase}/{year}/{sample}"
    print(outputDirectory)
    #####################################
    #Create Condor JDL file
    #####################################
    os.system("rm -f submit/MergeNtuple_{}_Job*.jdl".format(sample))
    os.system("rm -f log/MergeNtuple_{}_Job*".format(sample))


    jdl_file="submit/MergeNtuple_{}_{}.jdl".format( sample, maxjob)

    tmpCondorJDLFile = open(jdl_file,"w")

    tmpCondorJDLFile.write("Universe = vanilla \n")
    tmpCondorJDLFile.write("Executable = {} \n".format(executable))
    tmpCondorJDLFile.write("Arguments = {} {} $(ProcId) {} {} {} {}/ \n"\
                            .format(inputfilelist, filesPerJob, maxjob, outputDirectory, CMSSW_BASE, HOME))

    tmpCondorJDLFile.write("Log = log/MergeNtuple_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).log \n".format(sample, maxjob))
    tmpCondorJDLFile.write("Output = log/MergeNtuple_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).out \n".format(sample, maxjob))
    tmpCondorJDLFile.write("Error = log/MergeNtuple_{}_Job$(ProcId)_Of_{}_$(Cluster).$(Process).err \n".format(sample, maxjob))

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
