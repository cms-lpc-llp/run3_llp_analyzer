import os
import shutil
import glob

list_dirs = [
#    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_240310_CSCOnly_DNNs/V1p19/'
#    '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/Data2023/v10/*EXOCSCCluster*'
#    '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/Data2023/v10/Muon1-EXOCSCCluster_Run2023D-PromptReco-v1/',
#    '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/Data2023/v10/Muon1-EXOCSCCluster_Run2023C-PromptReco-v2/'

#    '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer23/v2/sixie/v2/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/'
    '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_240613_signals/V1p19_fixed/',
]

# make list of files 

for directory in list_dirs:
    for f in glob.glob(directory + '/**/*.root',recursive=True):
        os.system('echo ' + f + ' >> file_args.txt')
