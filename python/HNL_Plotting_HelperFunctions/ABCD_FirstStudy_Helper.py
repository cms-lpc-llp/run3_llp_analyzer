#!/usr/bin/env python3

import os

import coffea
import awkward as ak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot
import numpy as np
import hist
import yaml
import sys
import copy
import dask
from dask.distributed import Client
import gc

sys.path.append('.')



def return_clusterSize_dPhiClusterTau(events, passID = 'tauIsVLoose', failID = None, blind=True, hotspotCheck=False):
    '''
    Helper function for returning (computed) cluster size branch and dPhi(cluster, tau) branches with necessary preselections
    Cannot run unblinded without inverted TauID
    Baseline Preselection Applied:
    Pass Noise Filters
    Pass Trigger
    Exactly one reconstructed tau
    tau pass vLoose ID
    At least one reconstructed cluster (should happen 100% of the time because of trigger requirement)
    Reco cluster has >=100 hits, is in time,  and has no hits in ME11/12 (also should be near 100% of the time because trigger)
    Jet and Muon Veto pTs
    '''

    # if not invertTightID and not blind:
    #     print("Cannot unblind signal regiion if tau tight ID is not inverted")
    #     return
    
    #event level cuts
    noise_filters_cut = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
    trigger_cut = events.HLT_CscCluster100_PNetTauhPFJet10_Loose
    nTaus_cut = events.nTaus==1
    nClusters_cut = events.nCscRechitClusters>0

    events = events[(noise_filters_cut) & (trigger_cut) & (nTaus_cut) & (nClusters_cut)]

    # if invertTightID:
    #     tau_ID_cut = ak.all(((events.tauIsVLoose) & ~(events.tauIsTight)), axis=1)
    # else:
    #     tau_ID_cut = ak.all(((events.tauIsTight)), axis=1)

    if passID is not None and failID is not None:
        tau_ID_cut = ak.all(((events[passID]) & ~(events[failID])), axis=1)
    elif passID is not None:
        tau_ID_cut = ak.all(((events[passID])), axis=1)
    elif failID is not None:
        tau_ID_cut = ak.all((~(events[failID])), axis=1)
    else:
        tau_ID_cut = ak.ones_like((events.evtNum), axis=1, dtype=bool)
    
    events = events[tau_ID_cut]

    cluster_cuts = ((events.cscRechitClusterSize>=100) & (events.cscRechitClusterTimeWeighted > -5) & (events.cscRechitClusterTimeWeighted < 12.5) 
    & ((events.cscRechitClusterNRechitChamberMinus11 + events.cscRechitClusterNRechitChamberMinus12 + 
        events.cscRechitClusterNRechitChamberPlus11 + events.cscRechitClusterNRechitChamberPlus12)==0)
        & (events.cscRechitClusterMuonVetoPt<30) & (events.cscRechitClusterJetVetoPt<30))
    
    if blind:
        cluster_cuts = (cluster_cuts) & ~((events.cscRechitClusterSize>160) & (abs(events.cscRechitClusterPromptTauDeltaPhi)>2))
    
    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")

    if not hotspotCheck:
        clusterSize = events.cscRechitClusterSize[cluster_cuts]
        dPhi = abs(events.cscRechitClusterPromptTauDeltaPhi[cluster_cuts])
        return_tuple = ak.flatten(clusterSize[ak.num(clusterSize)>0]).compute(), ak.flatten(dPhi[ak.num(dPhi)>0]).compute() #remove empty entries

    else:
        clusterEta = events.cscRechitClusterEta[cluster_cuts]
        clusterPhi = events.cscRechitClusterPhi[cluster_cuts]
        return_tuple = ak.flatten(clusterEta[ak.num(clusterEta)>0]).compute(), ak.flatten(clusterPhi[ak.num(clusterPhi)>0]).compute() #remove empty entries

    
    client.close()

    return return_tuple

def return_clusterSize_dPhiClusterTau_allSelections(events, passID = 'tauIsVLoose', failID = None, blind=True, hotspotCheck=False):
    '''
    Helper function for returning (computed) cluster size branch and dPhi(cluster, tau) branches with necessary preselections
    Cannot run unblinded without inverted TauID
    Baseline Preselection Applied:
    Pass Noise Filters
    Pass Trigger
    Exactly one reconstructed tau
    tau pass vLoose ID
    At least one reconstructed cluster (should happen 100% of the time because of trigger requirement)
    Reco cluster has >=100 hits, is in time,  and has no hits in ME11/12 (also should be near 100% of the time because trigger)
    Jet and Muon Veto pTs

    ADDED MORE FOR THIS FUNCTION
    '''

    # if not invertTightID and not blind:
    #     print("Cannot unblind signal regiion if tau tight ID is not inverted")
    #     return
    
    #event level cuts
    noise_filters_cut = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
    trigger_cut = events.HLT_CscCluster100_PNetTauhPFJet10_Loose
    nTaus_cut = events.nTaus==1
    nClusters_cut = events.nCscRechitClusters>0
    MET_cut = events.Puppimet>30

    events = events[(noise_filters_cut) & (trigger_cut) & (nTaus_cut) & (nClusters_cut) & (MET_cut)]

    # if invertTightID:
    #     tau_ID_cut = ak.all(((events.tauIsVLoose) & ~(events.tauIsTight)) & (events.tauPt>20), axis=1)
    # else:
    #     tau_ID_cut = ak.all(((events.tauIsTight) & (events.tauPt>20)), axis=1)

    if passID is not None and failID is not None:
        tau_ID_cut = ak.all(((events[passID]) & ~(events[failID])), axis=1)
    elif passID is not None:
        tau_ID_cut = ak.all(((events[passID])), axis=1)
    elif failID is not None:
        tau_ID_cut = ak.all((~(events[failID])), axis=1)
    else:
        tau_ID_cut = ak.ones_like((events.evtNum), axis=1, dtype=bool)

    events = events[tau_ID_cut]

    cluster_cuts = ((events.cscRechitClusterSize>=100) & (events.cscRechitClusterTimeWeighted > -5) & (events.cscRechitClusterTimeWeighted < 12.5) 
    & ((events.cscRechitClusterNRechitChamberMinus11 + events.cscRechitClusterNRechitChamberMinus12 + 
        events.cscRechitClusterNRechitChamberPlus11 + events.cscRechitClusterNRechitChamberPlus12 
         + events.cscRechitClusterNRechitChamberPlus13 + events.cscRechitClusterNRechitChamberMinus13+
        events.cscRechitClusterNRechitChamberPlus21 + events.cscRechitClusterNRechitChamberPlus21 +
        events.cscRechitClusterNRechitChamberPlus22 + events.cscRechitClusterNRechitChamberMinus22
         )==0)  &
        (events.cscRechitClusterMuonVetoPt<30) & (events.cscRechitClusterJetVetoPt<30) & 
        (events.cscRechitClusterNStation10>1) &
        #(events.cscRechitClusterDNN_bkgMC_plusBeamHalo>=0.995) &
        (abs(events.cscRechitClusterPuppiMet_dPhi)<1.5) & (abs(events.cscRechitClusterPromptTauDeltaEta)<2))
    
    if blind:
        cluster_cuts = (cluster_cuts) & ~((events.cscRechitClusterSize>160) & (abs(events.cscRechitClusterPromptTauDeltaPhi)>2))
    
    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")

    #print(events.cscRechitClusterSize.compute())
    #print(cluster_cuts.compute())

    if not hotspotCheck:
        clusterSize = events.cscRechitClusterSize[cluster_cuts]
        dPhi = abs(events.cscRechitClusterPromptTauDeltaPhi[cluster_cuts])
        return_tuple = ak.flatten(clusterSize[ak.num(clusterSize)>0]).compute(), ak.flatten(dPhi[ak.num(dPhi)>0]).compute() #remove empty entries

    else:
        clusterEta = events.cscRechitClusterEta[cluster_cuts]
        clusterPhi =events.cscRechitClusterPhi[cluster_cuts]
        return_tuple = ak.flatten(clusterEta[ak.num(clusterEta)>0]).compute(), ak.flatten(clusterPhi[ak.num(clusterPhi)>0]).compute() #remove empty entries
    client.close()

    return return_tuple



def return_clusterSize_dPhiClusterTau_allSelectionsDNN(events, passID = 'tauIsVLoose', failID = None, DNN_cut=0.995, blind=True, hotspotCheck=False, additional_branches = []):
    '''
    Helper function for returning (computed) cluster size branch and dPhi(cluster, tau) branches with necessary preselections
    Cannot run unblinded without inverted TauID
    Baseline Preselection Applied:
    Pass Noise Filters
    Pass Trigger
    Exactly one reconstructed tau
    tau pass vLoose ID
    At least one reconstructed cluster (should happen 100% of the time because of trigger requirement)
    Reco cluster has >=100 hits, is in time,  and has no hits in ME11/12 (also should be near 100% of the time because trigger)
    Jet and Muon Veto pTs

    ADDED MORE FOR THIS FUNCTION
    '''

    # if not invertTightID and not blind:
    #     print("Cannot unblind signal regiion if tau tight ID is not inverted")
    #     return
    
    #event level cuts
    noise_filters_cut = (events.Flag_all) & (events.quot) & (events.Flag_ecalBadCalibFilter)
    trigger_cut = events.HLT_CscCluster100_PNetTauhPFJet10_Loose
    nTaus_cut = events.nTaus==1
    nClusters_cut = events.nCscRechitClusters>0
    MET_cut = events.Puppimet>30

    events = events[(noise_filters_cut) & (trigger_cut) & (nTaus_cut) & (nClusters_cut) & (MET_cut)]

    # if invertTightID:
    #     tau_ID_cut = ak.all(((events.tauIsVLoose) & ~(events.tauIsTight)) & (events.tauPt>20), axis=1)
    # else:
    #     tau_ID_cut = ak.all(((events.tauIsTight) & (events.tauPt>20)), axis=1)


    if passID is not None and failID is not None:
        tau_ID_cut = ak.all(((events[passID]) & ~(events[failID])), axis=1)
    elif passID is not None:
        tau_ID_cut = ak.all(((events[passID])), axis=1)
    elif failID is not None:
        tau_ID_cut = ak.all((~(events[failID])), axis=1)
    else:
        tau_ID_cut = ak.ones_like((events.evtNum), axis=1, dtype=bool)

    events = events[tau_ID_cut]

    cluster_cuts = ((events.cscRechitClusterSize>=100) & (events.cscRechitClusterTimeWeighted > -5) & (events.cscRechitClusterTimeWeighted < 12.5) 
    & ((events.cscRechitClusterNRechitChamberMinus11 + events.cscRechitClusterNRechitChamberMinus12 + 
        events.cscRechitClusterNRechitChamberPlus11 + events.cscRechitClusterNRechitChamberPlus12 
        #   + events.cscRechitClusterNRechitChamberPlus13 + events.cscRechitClusterNRechitChamberMinus13+
        #  events.cscRechitClusterNRechitChamberPlus21 + events.cscRechitClusterNRechitChamberPlus21 +
        #  events.cscRechitClusterNRechitChamberPlus22 + events.cscRechitClusterNRechitChamberMinus22
         )==0)  &
        (events.cscRechitClusterMuonVetoPt<30) & (events.cscRechitClusterJetVetoPt<30) & 
        #(events.cscRechitClusterNStation10>1) &
        (events.cscRechitClusterDNN_bkgMC_plusBeamHalo>=DNN_cut) &
        (abs(events.cscRechitClusterPuppiMet_dPhi)<1.5) & (abs(events.cscRechitClusterPromptTauDeltaEta)<2))
    
    if blind:
        cluster_cuts = (cluster_cuts) & ~((events.cscRechitClusterSize>160) & (abs(events.cscRechitClusterPromptTauDeltaPhi)>2))
    
    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")

    #print(events.cscRechitClusterSize.compute())
    #print(cluster_cuts.compute())

    if not hotspotCheck:
        clusterSize = events.cscRechitClusterSize[cluster_cuts]
        dPhi = abs(events.cscRechitClusterPromptTauDeltaPhi[cluster_cuts])
        return_tuple = ak.flatten(clusterSize[ak.num(clusterSize)>0]).compute(), ak.flatten(dPhi[ak.num(dPhi)>0]).compute() #remove empty entries

    else:
        clusterEta = events.cscRechitClusterEta[cluster_cuts]
        clusterPhi = events.cscRechitClusterPhi[cluster_cuts]
        return_tuple = ak.flatten(clusterEta[ak.num(clusterEta)>0]).compute(), ak.flatten(clusterPhi[ak.num(clusterPhi)>0]).compute() #remove empty entries
    

    for branch in additional_branches:
        if "cscCluster" in branch:
            new_branch = events[branch][cluster_cuts]
            return_tuple = return_tuple + (ak.flatten(new_branch[ak.num(new_branch)>0]).compute(),)

        elif "tau" in branch:
            new_branch = events[branch][ak.any(cluster_cuts, axis=1)]
            return_tuple = return_tuple + (ak.flatten(new_branch).compute(),)

        else: #event level
            new_branch = events[branch][ak.any(cluster_cuts, axis=1)]
            return_tuple = return_tuple + (new_branch.compute(),)
        

    client.close()

    return return_tuple
    

