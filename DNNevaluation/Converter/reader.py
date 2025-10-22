import ROOT
import argparse
from math import *
import os

# python reader.py --inputFile /eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p5_Removed2ClustersCut_CountZCSC/V1p19/MC_Summer22EE/v1/sixie/v5/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_Job0_of_371.root 

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Receive the parameters")
    parser.add_argument('--inputFile', action = 'store', type = str, dest = 'inputFile', help = 'Define the inputDir path')

    args = parser.parse_args()

#    variables = ['XXXRechitClusterSize', 'XXXRechitCluster_match_gLLP', 'XXXRechitClusterTime', 'XXXRechitClusterEta', 'XXXRechitClusterPhi', 'XXXRechitClusterNRechitChamberPlus11', 'XXXRechitClusterNRechitChamberMinus11', 'XXXRechitClusterNRechitChamberPlus12', 'XXXRechitClusterNRechitChamberMinus12', 'XXXRechitClusterXSpread', 'XXXRechitClusterYSpread', 'XXXRechitClusterZSpread', 'XXXRechitClusterXYSpread', 'XXXRechitClusterRSpread']

    #out_muons = open('/eos/cms/store/group/phys_muon/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_230412/txtFiles/' + args.inputFile.split('/')[-1].split('.')[0] + '.txt', 'w')
    out_muons = open('/eos/cms/store/group/phys_muon/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_240304_CSC/txtFiles_data/' + args.inputFile.split('/')[-1].split('.')[0] + '.txt', 'w')
    variables = ['XXXRechitClusterSize', 'XXXRechitCluster_match_gLLP', 'XXXRechitClusterTime', 'XXXRechitClusterEta', 'XXXRechitClusterPhi', 'XXXRechitClusterNRechitChamberPlus11', 'XXXRechitClusterNRechitChamberMinus11', 'XXXRechitClusterNRechitChamberPlus12', 'XXXRechitClusterNRechitChamberMinus12', 'XXXRechitClusterTimeSpread', 'XXXRechitClusterMuonVetoPt', 'XXXRechitClusterMuonVetoGlobal', 'XXXRechitClusterJetVetoLooseId', 'XXXRechitClusterJetVetoTightId', 'XXXRechitClusterJetVetoPt', 'XXXRechitClusternXY', 'XXXRechitClusternZ', 'XXXRechitClusterXSpread', 'XXXRechitClusterYSpread', 'XXXRechitClusterZSpread', 'XXXRechitClusterXYSpread', 'XXXRechitClusterRSpread', 'XXXRechitClusterMajorAxis', 'XXXRechitClusterMinorAxis', 'XXXRechitClusterSkewX', 'XXXRechitClusterSkewY', 'XXXRechitClusterSkewZ', 'XXXRechitClusterKurtX', 'XXXRechitClusterKurtY', 'XXXRechitClusterKurtZ']
#    variables = ['XXXRechitClusterSize', 'XXXRechitCluster_match_gLLP', 'XXXRechitClusterTime', 'XXXRechitClusterNRechitChamberPlus11', 'XXXRechitClusterNRechitChamberMinus11', 'XXXRechitClusterNRechitChamberPlus12', 'XXXRechitClusterNRechitChamberMinus12', 'XXXRechitClusterTimeSpread', 'XXXRechitClusterMuonVetoPt', 'XXXRechitClusterMuonVetoGlobal', 'XXXRechitClusterJetVetoPt', 'XXXRechitClusternXY', 'XXXRechitClusternZ', 'XXXRechitClusterXSpread', 'XXXRechitClusterYSpread', 'XXXRechitClusterZSpread', 'XXXRechitClusterXYSpread', 'XXXRechitClusterRSpread', 'XXXRechitClusterMajorAxis', 'XXXRechitClusterMinorAxis', 'XXXRechitClusterSkewX', 'XXXRechitClusterSkewY', 'XXXRechitClusterSkewZ', 'XXXRechitClusterKurtX', 'XXXRechitClusterKurtY', 'XXXRechitClusterKurtZ', 'XXXRechitClusterEtaPhiSpread', 'XXXRechitClusterEtaSpread', 'XXXRechitClusterPhiSpread', 'XXXRechitClusterDeltaRSpread']

#    out_muons = open('/eos/cms/store/group/phys_smp/fernanpe/displacedJetMuonAnalyzer_CA0p6_noMerging_231001/V1p19/txtFiles_gun/' + args.inputFile.split('/')[-1].split('.')[0] + '.txt', 'w')
#    out_muons = open('/eos/user/f/fernanpe/txtFiles_DifferentMassPoints/' + args.inputFile.split('/')[-1].split('.')[0] + '.txt', 'w')
#    out_muons = open('/afs/cern.ch/user/f/fernanpe/' + args.inputFile.split('/')[-1].split('.')[0] + '.txt', 'w')

    
    rootfile = ROOT.TFile.Open(args.inputFile, "READ")
    tree = rootfile.Get("MuonSystem")
            
    for event in tree:

        for i in range(eval('event.nCscRechitClusters')):
            #line = str(eval('event.runNum')) + '\t' + str(eval('event.lumiSec')) + '\t' + str(eval('event.evtNum')) +  '\t'
            #line += '1' + '\t'
            line = '1' + '\t'
            for var in variables:
                line += str(eval('event.' + var.replace('XXX', 'csc') + '[' + str(i) + ']')) + '\t'
            line += '-999' + '\t' # dtRechitCluster_match_RPCBx_dPhi0p5
            
            #nhits per station
            nhits_s1 = (float(eval('event.cscRechitClusterNRechitChamberPlus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus13' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus13' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            nhits_s2 = (float(eval('event.cscRechitClusterNRechitChamberPlus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus22' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus22' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            nhits_s3 = (float(eval('event.cscRechitClusterNRechitChamberPlus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus32' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus32' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            nhits_s4 = (float(eval('event.cscRechitClusterNRechitChamberPlus41' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus42' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus41' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus42' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            line += str(nhits_s1) + '\t' + str(nhits_s2) + '\t' + str(nhits_s3) + '\t' + str(nhits_s4) + '\t'

            #nhits per ring/wheel
            nhits_rw1 = (float(eval('event.cscRechitClusterNRechitChamberPlus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus41' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus41' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            nhits_rw2 = (float(eval('event.cscRechitClusterNRechitChamberPlus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus22' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus32' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus42' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus22' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus32' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus42' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            nhits_rw3 = (float(eval('event.cscRechitClusterNRechitChamberPlus13' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus13' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            line += str(nhits_rw1) + '\t' + str(nhits_rw2) + '\t' + str(nhits_rw3) + '\t' 





            out_muons.write(line + '\n')

        # for i in range(eval('event.nDtRechitClusters')):
        #     #line = str(eval('event.runNum')) + '\t' + str(eval('event.lumiSec')) + '\t' + str(eval('event.evtNum')) +  '\t'
        #     #line += '0' + '\t'
        #     line = '0' + '\t'
        #     for var in variables:
        #         if(('11' in var) or ('12' in var) or ('Time' in var)):
        #             line += '-999' + '\t' #ME11-12 (not in DTs, no timing)
        #         else:
        #             line += str(eval('event.' + var.replace('XXX', 'dt') + '[' + str(i) + ']')) + '\t'
        #     line += str(eval('event.dtRechitCluster_match_RPCBx_dPhi0p5' + '[' + str(i) + ']')) + '\t' # dtRechitCluster_match_RPCBx_dPhi0p5
            
        #     #nhits per station
        #     nhits_s1 = float(eval('event.dtRechitClusterNHitStation1' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
        #     nhits_s2 = float(eval('event.dtRechitClusterNHitStation2' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
        #     nhits_s3 = float(eval('event.dtRechitClusterNHitStation3' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
        #     nhits_s4 = float(eval('event.dtRechitClusterNHitStation4' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))

        #     line += str(nhits_s1) + '\t' + str(nhits_s2) + '\t' + str(nhits_s3) + '\t' + str(nhits_s4) + '\t'

        #     nhits_rw1 = float(eval('event.dtRechitClusterNHitWheel0' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
        #     nhits_rw2 = float(eval('event.dtRechitClusterNHitWheel1' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
        #     nhits_rw3 = float(eval('event.dtRechitClusterNHitWheel2' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))

        #     line += str(nhits_rw1) + '\t' + str(nhits_rw2) + '\t' + str(nhits_rw3) + '\t'

        #     out_muons.write(line + '\n')

    out_muons.close()
    rootfile.Close()
