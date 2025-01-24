#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
from collections import OrderedDict
from pathlib import Path


os.system("mkdir -p submit")
os.system("mkdir -p log")
executable = "MergeNtuple"
runExec = "runMergeNtuple_cache.sh"
filesPerJob = 1
list_path="/storage/af/user/christiw/login-1/christiw/LLP/Run3/CMSSW_10_6_30/src/run3_llp_analyzer/scripts_local/merge_ntuples_commands/"
#list_path="/storage/af/user/christiw/login-1/christiw/LLP/Run3/CMSSW_10_6_30/src/run3_llp_analyzer/scripts_local/merge_ntuples_commands/MC_Summer22/"

outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/MergedNtuples/Run3/V1p19/"
HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"
#executable = "llp_MuonSystem_CA_merge"
#list_path = "/storage/af/user/christiw/login-1/christiw/LLP/Run3/CMSSW_10_6_30/src/run3_llp_analyzer//lists/MergedNtuples/Run3/V1p19/"
#filesPerJob=10

datasetList = OrderedDict() #list of ntuples and point to the correct nano list file
samples = os.listdir(list_path)
#samples = Path(list_path).glob("**/*")

for sample_path in samples:
    #sample = str(sample_path.stem)
    sample = sample_path
    #if not "Muon0-EXOCSCCluster_Run2023C-PromptReco-v1" in sample:continue
    #if not "Cluster" in sample:continue

    if "list" in sample:continue
    if not "JetMET" in sample:continue
    if not "2024" in sample:continue
    print(sample)
    # if not "2022" in sample:continue
    #if not "ggH" in sample:continue
    print("Preparing workflow for dataset: " + sample + "\n")

    #inputfilelist  =  str(sample_path) 
    inputfilelist = list_path + '/' + sample
    # check if list files exist
    if not os.path.exists(inputfilelist):
        print("ntuple list file: " + inputfilelist + " does not exist. skipping.")
        continue
    print(inputfilelist) 
    cmd = ["awk", "END{print NR}",inputfilelist]
    nfiles = int(subprocess.check_output(cmd).decode("utf-8"))
    maxjob=int(nfiles/filesPerJob)+1
    mod=int(nfiles%filesPerJob)
    if mod == 0: maxjob -= 1

    year = sample[sample.find("Run")+3:sample.find("Run")+7]
    outputDirectory = f"{outputDirectoryBase}/{year}/{sample}"
    print(outputDirectory)
    job_list = []
    for i in range(maxjob):
        temp_file_name = list(Path("log/").glob(f"{executable}_{sample}_Job{i}_Of_{maxjob}*.err"))
        if len(temp_file_name)==1:temp_file_name = temp_file_name[0]
        else: #more than one files with the same name or couldn't find err file
            print("something is wrong: ",temp_file_name, f"{executable}_{sample}_Job{i}_Of_{maxjob}*.err")
            continue
        with open(temp_file_name, "r") as f1:
            last_line = f1.readlines()[-1]
        if (not "Unable to load sec.protocol plugin libXrdSecztn.so" in last_line) and \
        (not "Error in <TChain::SetBranchAddress>: unknown branch" in last_line) and (not "Error in <THashList::Delete>:" in last_line):job_list.append(i)
    print("Resubmitting: ", job_list)
    continue
    #####################################
    #Create Condor JDL file
    #####################################
    for job_number in job_list:

        os.system("rm -f submit/MergeNtuple_{}_Job{}*.jdl".format(sample, job_number))
        os.system("rm -f log/MergeNtuple_{}_Job{}*".format(sample, job_number))


        jdl_file="submit/MergeNtuple_{}_Job{}_of_{}.jdl".format( sample, job_number, maxjob)

        tmpCondorJDLFile = open(jdl_file,"w")

        tmpCondorJDLFile.write("Universe = vanilla \n")
        tmpCondorJDLFile.write("Executable = {} \n".format(runExec))
        tmpCondorJDLFile.write("Arguments = {} {} {} {} {} {} {}/ \n"\
                                .format(inputfilelist, filesPerJob, job_number, maxjob, outputDirectory, CMSSW_BASE, HOME))

        tmpCondorJDLFile.write("Log = log/MergeNtuple_{}_Job{}_Of_{}_$(Cluster).$(Process).log \n".format(sample, job_number, maxjob))
        tmpCondorJDLFile.write("Output = log/MergeNtuple_{}_Job{}_Of_{}_$(Cluster).$(Process).out \n".format(sample, job_number, maxjob))
        tmpCondorJDLFile.write("Error = log/MergeNtuple_{}_Job{}_Of_{}_$(Cluster).$(Process).err \n".format(sample, job_number, maxjob))

        tmpCondorJDLFile.write("+JobQueue=\"Normal\" \n")
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
        tmpCondorJDLFile.write("Queue 1 \n")
        tmpCondorJDLFile.close()

        os.system("condor_submit {} --batch-name {}".format(jdl_file, sample))
