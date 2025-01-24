#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys


datasets = [
## based on recommendation as of May 5 2024
#'/MET/Run2022A-22Sep2023-v1/NANOAOD',
#'/MET/Run2022B-22Sep2023-v1/NANOAOD',
#'/MET/Run2022C-22Sep2023-v1/NANOAOD',
#'/JetMET/Run2022C-22Sep2023-v1/NANOAOD',
#'/JetMET/Run2022D-22Sep2023-v1/NANOAOD',
#'/JetMET/Run2022E-22Sep2023-v1/NANOAOD',
#'/JetMET/Run2022F-22Sep2023-v2/NANOAOD',
#'/JetMET/Run2022G-22Sep2023-v2/NANOAOD',
#'/JetMET0/Run2023B-22Sep2023-v1/NANOAOD',
#'/JetMET0/Run2023C-22Sep2023_v1-v1/NANOAOD',
#'/JetMET0/Run2023C-22Sep2023_v2-v1/NANOAOD',
#'/JetMET0/Run2023C-22Sep2023_v3-v1/NANOAOD',
#'/JetMET0/Run2023C-22Sep2023_v4-v1/NANOAOD',
#'/JetMET0/Run2023D-22Sep2023_v1-v1/NANOAOD',
#'/JetMET0/Run2023D-22Sep2023_v2-v1/NANOAOD',
#'/JetMET1/Run2023B-22Sep2023-v2/NANOAOD',
#'/JetMET1/Run2023C-22Sep2023_v1-v1/NANOAOD',
#'/JetMET1/Run2023C-22Sep2023_v2-v1/NANOAOD',
#'/JetMET1/Run2023C-22Sep2023_v3-v1/NANOAOD',
#'/JetMET1/Run2023C-22Sep2023_v4-v1/NANOAOD',
#'/JetMET1/Run2023D-22Sep2023_v1-v1/NANOAOD',
#'/JetMET1/Run2023D-22Sep2023_v2-v1/NANOAOD',
#'/JetMET0/Run2024A-PromptReco-v1/NANOAOD',
#'/JetMET0/Run2024B-PromptReco-v1/NANOAOD',
#'/JetMET0/Run2024C-PromptReco-v1/NANOAOD',
#'/JetMET0/Run2024D-PromptReco-v1/NANOAOD',
#'/JetMET0/Run2024E-PromptReco-v1/NANOAOD',
#'/JetMET0/Run2024E-PromptReco-v2/NANOAOD',
#'/JetMET0/Run2024F-PromptReco-v1/NANOAOD',
#'/JetMET0/Run2024G-PromptReco-v1/NANOAOD',
#'/JetMET1/Run2024A-PromptReco-v1/NANOAOD',
#'/JetMET1/Run2024B-PromptReco-v1/NANOAOD',
#'/JetMET1/Run2024C-PromptReco-v1/NANOAOD',
#'/JetMET1/Run2024D-PromptReco-v1/NANOAOD',
#'/JetMET1/Run2024E-PromptReco-v1/NANOAOD',
#'/JetMET1/Run2024E-PromptReco-v2/NANOAOD',
#'/JetMET1/Run2024F-PromptReco-v1/NANOAOD',
'/JetMET0/Run2024G-PromptReco-v1/NANOAOD',
]


path = os.environ['CMSSW_BASE'] + '/src/run3_llp_analyzer/lists/nanoAOD/'
for name in datasets:
    
    if 'Run2022' in name: temp_path = path + '/2022/'
    elif 'Run2023' in name: temp_path = path + '/2023/'
    elif 'Run2024' in name: temp_path = path + '/2024/'
    if not os.path.exists(temp_path): os.makedirs(temp_path)
    outputFile = "{}/{}-{}.txt".format(temp_path, name.split('/')[1],name.split('/')[2])
    command = "dasgoclient -query=\"file dataset=" + name + " \" > temp.list"
    os.system(command)

    with open('temp.list', "r") as f:
        lines = f.readlines()
        for index, line in enumerate(lines):
            lines[index] = "root://cmsxrootd.fnal.gov/" + line.strip() + "\n"

    with open(outputFile, "w") as f:
        for line in lines:
            f.write(line)
    print("made list {}".format(outputFile))
