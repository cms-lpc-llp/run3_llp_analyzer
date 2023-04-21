from ROOT import TFile, TTree
import tensorflow as tf
from tensorflow import keras
import argparse
import numpy as np
import pandas as pd
from array import array

# python EvaluateDNN.py --in_file /eos/user/f/fernanpe/temp/ggH_HToSSTobbbb_MH-125_MS-15_CTau1000_13p6TeV_Job279_of_371.root
# python EvaluateDNN.py --in_file /eos/user/f/fernanpe/temp/DisplacedJet-EXOCSCCluster_Run2022E-PromptReco-v1_Job87_of_126.root

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Receive the parameters")
    parser.add_argument('--in_file', action = 'store', type = str, dest = 'in_file', help = 'input file')
    args = parser.parse_args()

    File = TFile.Open(args.in_file, 'update')
    tree = File.Get("MuonSystem")

    model = keras.models.load_model("training_CA0p6_NoMerging_FlatClusterSize0to500_230415.h5")

    nEntries = tree.GetEntries()

    dtRechitClusterDNN = array('f', [-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9])
    cscRechitClusterDNN = array('f', [-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9])
    branch1 = tree.Branch('dtRechitClusterDNN', dtRechitClusterDNN, 'dtRechitClusterDNN[30]/F')
    branch2 = tree.Branch('cscRechitClusterDNN', cscRechitClusterDNN, 'cscRechitClusterDNN[30]/F')


    for event in tree:
        
        # Evaluate the model for each cluster 

        for i in range(eval('event.nCscRechitClusters')):

            # Build DNN input variables

            frac_s1 = (float(eval('event.cscRechitClusterNRechitChamberPlus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus13' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus13' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            frac_s2 = (float(eval('event.cscRechitClusterNRechitChamberPlus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus22' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus22' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            frac_s3 = (float(eval('event.cscRechitClusterNRechitChamberPlus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus32' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus32' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            frac_s4 = (float(eval('event.cscRechitClusterNRechitChamberPlus41' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus42' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus41' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus42' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))

            frac_rw1 = (float(eval('event.cscRechitClusterNRechitChamberPlus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus41' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus11' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus21' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus31' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus41' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            frac_rw2 = (float(eval('event.cscRechitClusterNRechitChamberPlus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus22' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus32' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberPlus42' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus12' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus22' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus32' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus42' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            frac_rw3 = (float(eval('event.cscRechitClusterNRechitChamberPlus13' + '[' + str(i) + ']')) + float(eval('event.cscRechitClusterNRechitChamberMinus13' + '[' + str(i) + ']'))) / float(eval('event.cscRechitClusterSize' + '[' + str(i) + ']'))
            
            inputs = np.array([1, eval('event.cscRechitClusterXSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterYSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterZSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterXYSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterRSpread' + '[' + str(i) + ']'), frac_s1, frac_s2, frac_s3, frac_s4, frac_rw1, frac_rw2, frac_rw3], ndmin=2)
            inputs[np.isnan(inputs)] = -999.9 #replace nan by -9999.9 following the training criteria

            # Evaluate model

            cscRechitClusterDNN[i] = model.predict(inputs)[0]
            
        # Fill untouched array elements

        for i in range(eval('event.nCscRechitClusters'), len(cscRechitClusterDNN)):
            cscRechitClusterDNN[i] = -999.9


        # Repeat same procedure for CSC clusters

        for i in range(eval('event.nDtRechitClusters')):

            frac_s1 = float(eval('event.dtRechitClusterNHitStation1' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
            frac_s2 = float(eval('event.dtRechitClusterNHitStation2' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
            frac_s3 = float(eval('event.dtRechitClusterNHitStation3' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
            frac_s4 = float(eval('event.dtRechitClusterNHitStation4' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))

            frac_rw1 = float(eval('event.dtRechitClusterNHitWheel0' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
            frac_rw2 = float(eval('event.dtRechitClusterNHitWheel1' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))
            frac_rw3 = float(eval('event.dtRechitClusterNHitWheel2' + '[' + str(i) + ']')) / float(eval('event.dtRechitClusterSize' + '[' + str(i) + ']'))

            inputs = np.array([0, eval('event.dtRechitClusterXSpread' + '[' + str(i) + ']'), eval('event.dtRechitClusterYSpread' + '[' + str(i) + ']'), eval('event.dtRechitClusterZSpread' + '[' + str(i) + ']'), eval('event.dtRechitClusterXYSpread' + '[' + str(i) + ']'), eval('event.dtRechitClusterRSpread' + '[' + str(i) + ']'), frac_s1, frac_s2, frac_s3, frac_s4, frac_rw1, frac_rw2, frac_rw3], ndmin=2)
            inputs[np.isnan(inputs)] = -999.9 #replace nan by -9999.9 following the training criteria

            dtRechitClusterDNN[i] = model.predict(inputs)[0]

        for i in range(eval('event.nDtRechitClusters'), len(dtRechitClusterDNN)):
            dtRechitClusterDNN[i] = -999.9
        
        # Fill new branches

        branch1.Fill()
        branch2.Fill()


    File.Write("", TFile.kOverwrite)
    File.Close()
