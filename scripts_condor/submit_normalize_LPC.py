#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict

sample = sys.argv[1] # e.g. Muon0-Run2024B-PromptReco-v1
input_end = sys.argv[2] #e.g. 2024_Data (what to append to HNL_Tau_Search directory path for output)
log_name = sys.argv[3] # e.g. 2024_Data_Test (name of directory for condor submit and log files)



#os.system("mkdir -p submit")
#os.system("mkdir -p log")
executable = "scripts_condor/normalize_LPC.sh"

HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')

version = "v13"
#for year in ['2022', '2023']:
log_dir = "condor_job_output/log_{}".format(log_name)
submit_dir = "condor_job_output/submit_{}".format(log_name)
os.system("mkdir -p {}".format(log_dir))
os.system("mkdir -p {}".format(submit_dir))
inputDir = "/store/group/lpclonglived/amalbert/HNL_Tau_Search/{}/{}".format(input_end,sample)

os.system("eosls -f {}/ > fileList_normalization/inputfilelist_{}.txt".format(inputDir, sample))

outputDir=inputDir  + "/normalized/"
print("outputDir: ", outputDir)
os.system("eosmkdir {}".format(outputDir))

#### Run on all signal samples ####
## year/isData
datasetList = OrderedDict()
#samples = os.listdir(inputDir)
#for s in samples:
#if "normalize" in s: continue
#if 'Data' in inputDir: datasetList[s.replace('.txt', '')] = ["2018", "yes"]
#else: datasetList[s.replace('.txt', '')] = ["2018", "no"]
############
#os.system("eval `scram runtime -sh`")

print("input directory: " + inputDir)
#for sample in datasetList.keys():

print("Normalizing for dataset: " + sample + "\n")

#year = datasetList[sample][0]
isData="no"
if "Muon" in sample or "HNL" in sample:
    isData="yes"


#####################################
#Create Condor JDL file
#####################################
os.system("rm -f {}/normalize_{}.jdl".format(submit_dir,sample))
os.system("rm -rf {}*".format(log_dir))


jdl_file="{}/normalize_{}.jdl".format(submit_dir, sample)

tmpCondorJDLFile = open(jdl_file,"w")

tmpCondorJDLFile.write("Universe = vanilla \n")
tmpCondorJDLFile.write("Executable = {} \n".format(executable))
tmpCondorJDLFile.write("Arguments = {} {} {} {} {} {} {} \n".format(isData, sample, inputDir, outputDir, CMSSW_BASE, HOME, 1))

tmpCondorJDLFile.write("Log = {}/normalize_{}_$(Cluster).$(Process).log \n".format(log_dir,sample))
tmpCondorJDLFile.write("Output = {}/normalize_{}_$(Cluster).$(Process).out \n".format(log_dir,sample))
tmpCondorJDLFile.write("Error = {}/normalize_{}_$(Cluster).$(Process).err \n".format(log_dir,sample))


tmpCondorJDLFile.write("+JobQueue=\"Short\" \n")
tmpCondorJDLFile.write("RequestMemory = 4000 \n")
tmpCondorJDLFile.write("RequestCpus = 1 \n")
tmpCondorJDLFile.write("RequestDisk = 4 \n")

tmpCondorJDLFile.write("+RunAsOwner = True \n")
tmpCondorJDLFile.write("+InteractiveUser = true \n")
tmpCondorJDLFile.write("+SingularityImage = \"/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9\" \n")
tmpCondorJDLFile.write('+SingularityBindCVMFS = True \n')
tmpCondorJDLFile.write("run_as_owner = True \n")
#tmpCondorJDLFile.write("x509userproxy = {}/x509_proxy \n".format(HOME))
tmpCondorJDLFile.write("should_transfer_files = YES \n")
tmpCondorJDLFile.write("when_to_transfer_output = ON_EXIT \n")
tmpCondorJDLFile.write("Transfer_Input_Files = fileList_normalization/inputfilelist_{}.txt, data/xSections.dat, normalization_input.txt, NormalizeNtuple, Tarballs_for_LPC_Condor/SSLTarball.tar.gz\n".format(sample))
tmpCondorJDLFile.write("Queue 1 \n")
tmpCondorJDLFile.close()

os.system("condor_submit {} --batch-name {}".format(jdl_file, "normalize_" + sample))
