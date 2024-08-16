from ROOT import TFile, TTree
import tensorflow as tf
from tensorflow import keras
import argparse
import numpy as np
import pandas as pd
from array import array
import sys

# python EvaluateDNN.py --in_file /eos/user/f/fernanpe/temp/DYto2Mu_MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8_Job9_of_70.root

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = "Receive the parameters")
    parser.add_argument('--in_file', action = 'store', type = str, dest = 'in_file', help = 'input file')
    args = parser.parse_args()

    File = TFile.Open(args.in_file, 'update')
    tree = File.Get("MuonSystem")
    
    # Already filled!
    if 'cscRechitClusterDNN_bkgMC_plusBeamHalo' in [key.GetName() for key in tree.GetListOfBranches()]:
        sys.exit(0)

    #model_bkgMC = keras.models.load_model("training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_240505.h5")
    model_bkgMC_plusBeamHalo = keras.models.load_model("training_CA0p6_NoMerging_WeightedClusterSize_bkgMC_CSCOnly_adversarial_PlusBeamHalo_240510.h5")
    #model_bkgOOTData = keras.models.load_model("training_CA0p6_NoMerging_WeightedClusterSize_bkgDATA_CSCOnly_240510.h5")

    nEntries = tree.GetEntries()

    cscRechitClusterDNN_bkgMC = array('f', [-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9])
    cscRechitClusterDNN_bkgMC_plusBeamHalo = array('f', [-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9])
    cscRechitClusterDNN_bkgOOTData = array('f', [-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9, -999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9,-999.9])

    #branch1 = tree.Branch('cscRechitClusterDNN_bkgMC', cscRechitClusterDNN_bkgMC, 'cscRechitClusterDNN_bkgMC[30]/F')
    branch2 = tree.Branch('cscRechitClusterDNN_bkgMC_plusBeamHalo', cscRechitClusterDNN_bkgMC_plusBeamHalo, 'cscRechitClusterDNN_bkgMC_plusBeamHalo[30]/F')
    #branch3 = tree.Branch('cscRechitClusterDNN_bkgOOTData', cscRechitClusterDNN_bkgOOTData, 'cscRechitClusterDNN_bkgOOTData[30]/F')


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
            
            inputs = np.array([eval('event.cscRechitClusterXSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterYSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterZSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterXYSpread' + '[' + str(i) + ']'), eval('event.cscRechitClusterRSpread' + '[' + str(i) + ']'),  eval('event.cscRechitClusterSkewX' + '[' + str(i) + ']'),  eval('event.cscRechitClusterSkewY' + '[' + str(i) + ']'), frac_s1, frac_s2, frac_s3, frac_s4, frac_rw1, frac_rw2, frac_rw3], ndmin=2)
            inputs[np.isnan(inputs)] = -999.9 #replace nan by -9999.9 following the training criteria

            # Evaluate model

            #cscRechitClusterDNN_bkgMC[i] = model_bkgMC.predict(inputs)[0]
            cscRechitClusterDNN_bkgMC_plusBeamHalo[i] = model_bkgMC_plusBeamHalo.predict(inputs)[0]
            #cscRechitClusterDNN_bkgOOTData[i] = model_bkgOOTData.predict(inputs)[0]

        # Fill untouched array elements
        for i in range(eval('event.nCscRechitClusters'), len(cscRechitClusterDNN_bkgMC)):
            #cscRechitClusterDNN_bkgMC[i] = -999.9
            cscRechitClusterDNN_bkgMC_plusBeamHalo[i] = -999.9
            #cscRechitClusterDNN_bkgOOTData[i] = -999.9
        
        # Fill new branches

        #branch1.Fill()
        branch2.Fill()
        #branch3.Fill()


    File.Write("", TFile.kOverwrite)
    File.Close()
