#!/usr/bin/python

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
decay_list = ['4B','4D','4Tau']
mass = [1, 7, 15, 23, 30, 40, 55] #GeV
ct_list = ['0p1','1','10','100','1000','10000','100000',]
base_path = os.environ['CMSSW_BASE'] + '/src/run3_llp_analyzer/lists/nanoAOD/'
for p in tag.keys():
    for decay in decay_list:
        for m in mass:
            for ct in ct_list:
                if decay == '4B' and m < 15:continue
                if decay == '4D' and m < 7:continue
                list_path = '{}/MC_{}/'.format(base_path, p)
                print(list_path)
                if not os.path.exists(list_path): os.makedirs(list_path)
                sample_name = "ggH_Hto2Sto{}_MH-125-MS-{}-ctauS-{}_TuneCP5_13p6TeV_powheg-pythia8".format(decay, m, ct)
                outputFile = list_path + sample_name + ".txt"
                sample_name_das = "/{}/{}/NANOAODSIM".format(sample_name, tag[p])
                print(sample_name_das)

                command = "dasgoclient -query=\"dataset=" + sample_name_das + "\""
                proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
                (dataset, err) = proc.communicate()
                dataset = dataset.decode("utf-8")
                dataset=str(dataset).replace('\n','')

                command_getfiles = "dasgoclient -query=\"file dataset=" + dataset + " \" > temp.list"
                print(command_getfiles)
                dataset=os.system(command_getfiles)
                print(dataset)
                with open('temp.list', "r") as f:
                    lines = f.readlines()
                    print("Number of files: ", len(lines))
                    for index, line in enumerate(lines):
                        lines[index] = "root://cmsxrootd.fnal.gov/" + line.strip() + "\n"
        
                with open(outputFile, "w") as f:
                    for line in lines:
                        f.write(line)
                print("made list {}".format(outputFile))
