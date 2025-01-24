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
executable = "normalize.sh"

HOME = os.getenv('HOME')
CMSSW_BASE = os.getenv('CMSSW_BASE')

version = "v19"
for year in ['2024']:
    inputDir = f"/storage/af/group/phys_exotica/delayedjets/displacedJetMuonAnalyzer/Run3/V1p19/Data{year}/{version}/"
    outputDir=inputDir  + "/normalized/"
    
    
    #### Run on all signal samples ####
    ## year/isData
    datasetList = OrderedDict()
    samples = os.listdir(inputDir)
    for s in samples:
        if "normalize" in s: continue
        if 'Data' in inputDir: datasetList[s.replace('.txt', '')] = ["2018", "yes"]
        else: datasetList[s.replace('.txt', '')] = ["2018", "no"]
    ############
    os.system("eval `scram runtime -sh`")
    
    print("input directory: " + inputDir)
    for sample in datasetList.keys():
    
        print("Normalizing for dataset: " + sample + "\n")
    
        year = datasetList[sample][0]
        isData = datasetList[sample][1]
    
        #####################################
        #Create Condor JDL file
        #####################################
        os.system("rm -f submit/normalize_{}*.jdl".format(sample))
        os.system("rm -f log/normalize_{}*".format(sample))
    
    
        jdl_file="submit/normalize_{}.jdl".format(sample)
    
        tmpCondorJDLFile = open(jdl_file,"w")
    
        tmpCondorJDLFile.write("Universe = vanilla \n")
        tmpCondorJDLFile.write("Executable = {} \n".format(executable))
        tmpCondorJDLFile.write("Arguments = {} {} {} {} {} {} {} \n".format(isData, sample, inputDir, outputDir, CMSSW_BASE, HOME, 1))
    
        tmpCondorJDLFile.write("Log = log/normalize_{}_$(Cluster).$(Process).log \n".format(sample))
        tmpCondorJDLFile.write("Output = log/normalize_{}_$(Cluster).$(Process).out \n".format(sample))
        tmpCondorJDLFile.write("Error = log/normalize_{}_$(Cluster).$(Process).err \n".format(sample))
    
    
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
        tmpCondorJDLFile.write("Queue 1 \n")
        tmpCondorJDLFile.close()
    
        os.system("condor_submit {} --batch-name {}".format(jdl_file, "normalize_" + sample))
