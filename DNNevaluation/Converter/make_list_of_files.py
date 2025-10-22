import os
import shutil
import glob

# list_dirs = [
# #    '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer22EE/v6/fernanpe/v10/'
#     #'/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v5/fernanpe/v10/',
# #    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v4/fernanpe/v10/QCD_PT-1to20_noPU_Filter_nMuonClusters_gt0/',
#     #'/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v3/fernanpe/v10/QCD_PT-20to40_Filter_nMuonClusters_gt0/'
# #    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v2/fernanpe/v10/QCD_PT-0p5to20_Filter_nMuonClusters_gt0/'
# #    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22/v1/sixie/v10/'
# ]

list_dirs = [
    #'/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer23/v1/sixie/v10/', #DY
    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/Data2022/v10/DisplacedJet-EXOCSCCluster_Run2022E-PromptReco-v1/', #data
    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/Data2022/v10/DisplacedJet-EXOCSCCluster_Run2022F-PromptReco-v1/', #data
    '/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/Data2022/v10/DisplacedJet-EXOCSCCluster_Run2022G-PromptReco-v1/', #data
    #'/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v1/sixie/v10/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/', #signal
    #'/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v3/fernanpe/v10/QCD_PT-20to40_Filter_nMuonClusters_gt0/', #for PU
    #'/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22EE/v2/fernanpe/v10/QCD_PT-0p5to20_Filter_nMuonClusters_gt0/', #for PU
    #'/eos/user/f/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/V1p19/MC_Summer22/v1/sixie/v10/', #different mass points
    # '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer22EE/v6/fernanpe/v10/Gun_ckaon_PT-1to20_negOE_noPU_nMuonClusters_gt0/',
    # '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer22EE/v6/fernanpe/v10/Gun_cpion_PT-1to20_negOE_noPU_nMuonClusters_gt0/',
    # '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer22EE/v6/fernanpe/v10/Gun_cpion_PT-1to20_posOE_noPU_nMuonClusters_gt0/',
    # '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer22EE/v6/fernanpe/v10/Gun_k0L_PT-1to20_negOE_noPU_nMuonClusters_gt0/',
    # '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer22EE/v6/fernanpe/v10/Gun_k0L_PT-1to20_posOE_noPU_nMuonClusters_gt0/',
    # '/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/MC_Summer22EE/v6/fernanpe/v10/Gun_ckaon_PT-1to20_posOE_noPU_nMuonClusters_gt0/',

]

# make list of files 

for directory in list_dirs:
    for f in glob.glob(directory + '/**/*.root',recursive=True):
        os.system('echo ' + f + ' >> file_args.txt')
