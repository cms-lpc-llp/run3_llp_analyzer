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
executable = "runMergeNtuple.sh"
NtuplePerJob = 1
NanoPerJob = 1
year = '2022'
ntupler_version = 'V1p19/Data{}/'.format(year)



analyzer_version = 'v1'
outputDirectoryBase="/storage/af/group/phys_exotica/delayedjets/MergedNtuples/Run3/V1p19/{}/{}/".format(ntupler_version, analyzer_version)
HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')
Analyzer_DIR = CMSSW_BASE+"/src/run3_llp_analyzer/"

ntupleListDir = Analyzer_DIR + "lists/displacedJetMuonNtuple/{}/".format(ntupler_version)
nanoAODListDir = Analyzer_DIR + "lists/nanoAOD/{}/".format(year)

datasetList = OrderedDict() #list of ntuples and point to the correct nano list file
samples = os.listdir(ntupleListDir)


for s in samples:
    if not "JetMET" in s:continue # run over JETMET dataset
    era = s[s.find('Run'):s.find('Run')+8]
    for nano in os.listdir(nanoAODListDir):
        if era in s:
            datasetList[s.replace('.txt', '')] = nano.replace('.txt', '')
            break

############
for k,v in datasetList.items():print(k,v)

for sample in datasetList.keys():

    print("Preparing workflow for dataset: " + sample + "\n")

    inputNtupleList  = ntupleListDir + sample +'.txt'
    inputNanoList  = nanoAODListDir + datasetList[sample] +'.txt'

    # check if list files exist
    if not os.path.exists(inputNtupleList):
        print("ntuple list file: " + inputNtupleList + " does not exist. skipping.")
        continue
    if not os.path.exists(inputNanoList):
        print("nano list file: " + inputNanoList + " does not exist. skipping.")
        continue

    # split ntuple files based on NtuplePerJob
    cmd = ["awk", "END{print NR}",inputNtupleList]
    nNtupleFile = int(subprocess.check_output(cmd).decode("utf-8"))
    maxNtupleJob=int(nNtupleFile/NtuplePerJob)+1
    mod=int(nNtupleFile%NtuplePerJob)
    if mod == 0: maxNtupleJob -= 1

    # split NANO files based on NanoPerJob
    cmd = ["awk", "END{print NR}",inputNanoList]
    nNanoFile = int(subprocess.check_output(cmd).decode("utf-8"))
    maxNanoJob=int(nNanoFile/NanoPerJob)+1
    mod=int(nNanoFile%NanoPerJob)
    if mod == 0: maxNanoJob -= 1


    total_job = maxNtupleJob * maxNanoJob
    print("Number of Ntuple jobs: {}".format(maxNtupleJob))
    print("Number of nanoAOD jobs: {}".format(maxNanoJob))
    print("Total number of jobs: {}".format(total_job))

    outputDirectory = outputDirectoryBase + sample + "/"

    #####################################
    #Create Condor JDL file
    #####################################
    os.system("rm -f submit/MergeNtuple_{}_*.jdl".format(sample))
    os.system("rm -f log/MergeNtuple_{}_*".format(sample))


    jdl_file="submit/MergeNtuple_{}_{}.jdl".format(sample, total_job)

    tmpCondorJDLFile = open(jdl_file,"w")

    tmpCondorJDLFile.write("Universe = vanilla \n")
    tmpCondorJDLFile.write("Executable = {} \n".format(executable))
    

    # tmpCondorJDLFile.write("maxNtupleJob = {} \n".format(maxNtupleJob))
    tmpCondorJDLFile.write("maxNanoJob = {} \n".format(maxNanoJob))
    tmpCondorJDLFile.write("ntuple_index = ($(ProcId) / $(maxNanoJob)) \n")
    tmpCondorJDLFile.write("nano_index = ($(ProcId) % $(maxNanoJob)) \n")
    tmpCondorJDLFile.write("Arguments = {} {} $INT(ntuple_index) {} {} {} $INT(nano_index) {} $(ProcId) {} {} {}/ \n"\
                            .format(inputNtupleList, NtuplePerJob, maxNtupleJob, inputNanoList,  NanoPerJob, maxNanoJob, outputDirectory, CMSSW_BASE, HOME))
    


    log_string = "MergeNtuple_{}_Job$(ProcId)_Of_{}_NtupleJob$INT(ntuple_index)_Of_{}_NanoJob$INT(nano_index)_Of_{}_$(Cluster).$(Process)".format(sample, total_job, maxNtupleJob, maxNanoJob)
    tmpCondorJDLFile.write("Log = log/{}.log \n".format(log_string))
    tmpCondorJDLFile.write("Output = log/{}.out \n".format(log_string))
    tmpCondorJDLFile.write("Error = log/{}.err \n".format(log_string))

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
    tmpCondorJDLFile.write("Queue {} \n".format(total_job))
    tmpCondorJDLFile.close()

    os.system("condor_submit {} --batch-name {}".format(jdl_file, sample))
