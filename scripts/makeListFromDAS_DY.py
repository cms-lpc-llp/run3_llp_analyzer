#!/usr/bin/python
import subprocess
import os
import datetime
import time
import subprocess
import numpy as np
import glob
import sys
tag = {'Summer22EE':'Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_*',\
        'Summer22': 'Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_*',\
        'Summer23BPix':'Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_*',\
        'Summer23': 'Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_*'}
sample_list = [
"DYto2Mu_MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-1500to2500_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-2500to4000_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-4000to6000_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-6000_TuneCP5_13p6TeV_powheg-pythia8",
"DYto2Mu_MLL-800to1500_TuneCP5_13p6TeV_powheg-pythia8",
]

base_path = os.environ['CMSSW_BASE'] + '/src/run3_llp_analyzer/lists/nanoAOD/'
for p in tag.keys():
    for sample_name in sample_list:
        list_path = '{}/MC_{}/'.format(base_path, p)
        print(list_path)
        if not os.path.exists(list_path): os.makedirs(list_path)
        outputFile = list_path + sample_name + ".txt"
        sample_name_das = "/{}/{}/NANOAODSIM".format(sample_name, tag[p])
        print(sample_name_das) 
        command = "dasgoclient -query=\"dataset=" + sample_name_das + "\"" 
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (dataset, err) = proc.communicate()
#        print("program output:", out)

        dataset = dataset.decode("utf-8")
        dataset=str(dataset).replace('\n','')
        command_getfiles = "dasgoclient -query=\"file dataset=" + dataset  + " \" > temp.list"
        os.system(command_getfiles)
        print(command_getfiles)
        with open('temp.list', "r") as f:
            lines = f.readlines()
            print("Number of files: ", len(lines))
            for index, line in enumerate(lines):
                lines[index] = "root://cmsxrootd.fnal.gov/" + line.strip() + "\n"
    
        with open(outputFile, "w") as f:
            for line in lines:
                f.write(line)
        print("made list {}".format(outputFile))
