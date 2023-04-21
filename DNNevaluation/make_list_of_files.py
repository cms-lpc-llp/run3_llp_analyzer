import os
import shutil
import glob

list_dirs = [
    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/Data2022/v10/DisplacedJet-EXOCSCCluster_Run2022E-PromptReco-v1/',
    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/Data2022/v10/DisplacedJet-EXOCSCCluster_Run2022F-PromptReco-v1/',
    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/Data2022/v10/DisplacedJet-EXOCSCCluster_Run2022G-PromptReco-v1/',
    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v1/sixie/v10/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/',
]

# make list of files 

for directory in list_dirs:
    for f in glob.glob(directory + '/**/*.root',recursive=True):
        os.system('echo ' + f + ' >> file_args.txt')
