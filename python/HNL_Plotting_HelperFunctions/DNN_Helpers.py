#!/usr/bin/env python3

import sys
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
import pandas as pd
import dask_awkward as dak
import dask
import dask.dataframe as dd
from dask.distributed import Client
import modeling_cut_lookup


def get_DNN_array(events, clusterSizeCut=False, type='tau', isData=True):
    '''
    Function to return just DNN array from files or list of files in data or MC
    Cuts to apply (after noise filters):
    1. Pass Trigger
    2. At least one reconstructed cluster with >= 100 hits
    3. Cluster is in time
    4. Cluster has no hits in ME11/12
    5. Cluster Passes Muon and Jet Vetos
    5. (Optional) cluster size is larger than 160
    '''
    if isData:
        events = abs(events.cscRechitClusterPromptTauDeltaEta)>2
    if type=='tau':
        trigger_path = "HLT_CscCluster100_PNetTauhPFJet10_Loose"
        muonVetoPtThreshold = 30
    elif type=='e':
        trigger_path = "HLT_CscCluster100_Ele5"
        muonVetoPtThreshold = 30
    elif type=='mu':
        trigger_path = "HLT_CscCluster100_Mu5"
        muonVetoPtThreshold = 0.0000001
    else:
        raise ValueError("Invalid type specified. Choose from 'tau', 'e', or 'mu'.")


    noise_filters_cut = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
    trigger_cut = events[trigger_path]
    nClusters_cut = events.nCscRechitClusters>0

    events = events[(noise_filters_cut) & (trigger_cut) & (nClusters_cut)]

    cluster_cuts = ((events.cscRechitClusterSize>=100) & (events.cscRechitClusterTimeWeighted > -5) & (events.cscRechitClusterTimeWeighted < 12.5) 
    & ((events.cscRechitClusterNRechitME1112)==0)
        & (events.cscRechitClusterMuonVetoPt<muonVetoPtThreshold) & (events.cscRechitClusterJetVetoPt<10))

    if clusterSizeCut:
        cluster_cuts = (cluster_cuts) & (events.cscRechitClusterSize>=160)

    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, 
                local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
    
    DNN_Scores = events.cscRechitClusterDNN_bkgMC_plusBeamHalo[cluster_cuts]
    DNN_Scores = ak.flatten(DNN_Scores[ak.num(DNN_Scores)>0]).compute()
    
    client.close()

    return DNN_Scores


def signalEff_bkgEff(DNN_Scores_Signal, DNN_Scores_Data, sampling_scores):
    '''
    From two (computed) DNN awkward arrays, return list of signal efficiencies,
    data efficiencies, signal eff/bkg eff for the specified sampling scores
    '''
    signal_len = len(DNN_Scores_Signal)
    data_len = len(DNN_Scores_Data)
    effs_signal, effs_data, eff_ratio = [],[],[]
    for score in sampling_scores:
        eff_signal = ak.count_nonzero(DNN_Scores_Signal>=score)/signal_len
        eff_data = ak.count_nonzero(DNN_Scores_Data>=score)/data_len
        effs_signal.append(eff_signal)
        effs_data.append(eff_data)
        eff_ratio.append(eff_signal/eff_data)
    return effs_signal, effs_data, eff_ratio


def get_stationVeto_eff(events, clusterSizeCut=False):
    '''
    Function to return efficiency of station veto + NStation10>1 selections
    Initial Preselections (after noise filters):
    1. Pass Trigger
    2. At least one reconstructed cluster with >= 100 hits
    3. Cluster is in time
    4. Cluster has no hits in ME11/12
    5. Cluster Passes Muon and Jet Vetos
    5. (Optional) cluster size is larger than 160
    '''

    noise_filters_cut = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
    trigger_cut = events.HLT_CscCluster100_PNetTauhPFJet10_Loose
    nClusters_cut = events.nCscRechitClusters>0

    events = events[(noise_filters_cut) & (trigger_cut) & (nClusters_cut)]

    cluster_cuts = ((events.cscRechitClusterSize>=100) & (events.cscRechitClusterTimeWeighted > -5) & (events.cscRechitClusterTimeWeighted < 12.5) 
    & ((events.cscRechitClusterNRechitChamberMinus11 + events.cscRechitClusterNRechitChamberMinus12 + 
        events.cscRechitClusterNRechitChamberPlus11 + events.cscRechitClusterNRechitChamberPlus12)==0)
        & (events.cscRechitClusterMuonVetoPt<30) & (events.cscRechitClusterJetVetoPt<30))

    if clusterSizeCut:
        cluster_cuts = (cluster_cuts) & (events.cscRechitClusterSize>=160)

    #denom = ak.count_nonzero(cluster_cuts)
    
    num_cuts = (cluster_cuts) & ((events.cscRechitClusterNRechitChamberMinus13 + events.cscRechitClusterNRechitChamberPlus13 + 
        events.cscRechitClusterNRechitChamberPlus21 + events.cscRechitClusterNRechitChamberMinus21 +
        events.cscRechitClusterNRechitChamberPlus22 + events.cscRechitClusterNRechitChamberMinus22)==0) & (events.cscRechitClusterNStation10>1)

    #num = ak.count_nonzero(num_cuts)

    #print(cluster_cuts.compute())
    #print(num_cuts.compute())

    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, 
                local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
    
    eff = (ak.count_nonzero(ak.flatten(num_cuts))/ak.count_nonzero(ak.flatten(cluster_cuts))).compute()
    
    client.close()

    return eff


def get_DNN_ClusterSize_arrays_Tony(events, type='tau', isData=True):
    '''
    Function to return DNN and cluster size arrays from files or list of files in data or MC
    Invert the dPhi(cluster, MET) to stay in validation region
    Cuts to apply (after noise filters):
    1. Pass Trigger
    2. At least one reconstructed cluster with >= 100 hits
    3. Cluster is in time
    4. Cluster timespread < 25
    5. Cluster has no hits in ME11/12
    6. Cluster Passes Muon and Jet Vetos from Double cluster
   
    '''
    
    if type=='tau':
        trigger_path = "HLT_CscCluster100_PNetTauhPFJet10_Loose"
        muonVetoPtThreshold = 30
    elif type=='e':
        trigger_path = "HLT_CscCluster100_Ele5"
        muonVetoPtThreshold = 30
    elif type=='mu':
        trigger_path = "HLT_CscCluster100_Mu5"
        muonVetoPtThreshold = 0.0000001
    else:
        raise ValueError("Invalid type specified. Choose from 'tau', 'e', or 'mu'.")

    #trigger_path = "HLT_CscCluster100_Nominal"
    muonVetoPtThreshold = 30
    noise_filters_cut = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
    trigger_cut = events[trigger_path]
    nClusters_cut = events.nCscRechitClusters==1

    events = events[(noise_filters_cut) & (trigger_cut) & (nClusters_cut)]

    cluster_cuts = ((events.cscRechitClusterSize>=100) & (events.cscRechitClusterSize<=400)&(events.cscRechitClusterTimeWeighted > -5) & (events.cscRechitClusterTimeWeighted < 12.5) 
    & ((events.cscRechitClusterNRechitME1112)==0) & (events.cscRechitClusterTimeSpreadWeightedAll<25) &
    (events.cscRechitClusterMuonVetoPt<muonVetoPtThreshold) & ~(events.cscRechitClusterMuonVetoLooseId) &(events.cscRechitClusterJetVetoPt<30))

    if not isData:
        cluster_cuts = (cluster_cuts) & (events.cscRechitCluster_match_gLLP)
    else:
        cluster_cuts = (cluster_cuts) & (abs(events.cscRechitClusterMet_dPhi)>2)
    # if clusterSizeCut:
    #     cluster_cuts = (cluster_cuts) & (events.cscRechitClusterSize>=160)

    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, 
                local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
    
    DNN_Scores = events.cscRechitClusterDNN_bkgMC_plusBeamHalo[cluster_cuts]
    DNN_Scores = ak.flatten(DNN_Scores[ak.num(DNN_Scores)>0]).compute()

    cluster_size = events.cscRechitClusterSize[cluster_cuts]
    cluster_size = ak.flatten(cluster_size[ak.num(cluster_size)>0]).compute()
    
    client.close()

    return DNN_Scores, cluster_size

def signalEff_bkgEff_Tony(DNN_Scores_Signal, DNN_Scores_Data, ClusterSize_Signal, ClusterSize_Data,sampling_scores, cluster_size_cut):
    '''
    From two (computed) DNN awkward arrays and cluster size arrays, return list of signal efficiencies,
    data efficiencies, signal eff/bkg eff for the specified sampling scores
    '''
    signal_len = len(DNN_Scores_Signal)
    data_len = len(DNN_Scores_Data)
    effs_signal, effs_data, eff_ratio = [],[],[]
    score_MC_clustersize = modeling_cut_lookup._lookup_mc_cut(cluster_size_cut, "cluster_size") [1]
    print(score_MC_clustersize)
    for score in sampling_scores:
        score_MC_DNN = modeling_cut_lookup._lookup_mc_cut(score, "dnn")[1]
        if score<0.96:score_MC_DNN=score
        #print(score_MC_DNN)
        eff_signal = ak.count_nonzero((DNN_Scores_Signal>=score_MC_DNN)&(ClusterSize_Signal>=score_MC_clustersize))/signal_len
        eff_data = ak.count_nonzero((DNN_Scores_Data>=score)&(ClusterSize_Data>=cluster_size_cut))/data_len
        effs_signal.append(eff_signal)
        effs_data.append(eff_data)
        eff_ratio.append(eff_signal/eff_data)
    return effs_signal, effs_data, eff_ratio


def get_DNN_ClusterSize_ClustereZ_arrays(events, type='tau', isData=True):
    '''
    Function to return DNN and cluster size arrays from files or list of files in data or MC
    Invert the dPhi(cluster, MET) to stay in validation region
    Cuts to apply (after noise filters):
    1. Pass Trigger
    2. At least one reconstructed cluster with >= 100 hits
    3. Cluster is in time
    4. Cluster timespread < 25
    5. Cluster has no hits in ME11/12
    6. Cluster Passes Muon and Jet Vetos from Double cluster
   
    '''
    
    if type=='tau':
        trigger_path = "HLT_CscCluster100_PNetTauhPFJet10_Loose"
        muonVetoPtThreshold = 30
    elif type=='e':
        trigger_path = "HLT_CscCluster100_Ele5"
        muonVetoPtThreshold = 30
    elif type=='mu':
        trigger_path = "HLT_CscCluster100_Mu5"
        muonVetoPtThreshold = 0.0000001
    else:
        raise ValueError("Invalid type specified. Choose from 'tau', 'e', or 'mu'.")

    #trigger_path = "HLT_CscCluster100_Nominal"
    muonVetoPtThreshold = 30
    noise_filters_cut = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
    trigger_cut = events[trigger_path]
    nClusters_cut = events.nCscRechitClusters==1

    events = events[(noise_filters_cut) & (trigger_cut) & (nClusters_cut)]

    cluster_cuts = ((events.cscRechitClusterSize>=100) & (events.cscRechitClusterSize<=400)&(events.cscRechitClusterTimeWeighted > -5) & (events.cscRechitClusterTimeWeighted < 12.5) 
    & ((events.cscRechitClusterNRechitME1112)==0) & (events.cscRechitClusterTimeSpreadWeightedAll<25) &
    (events.cscRechitClusterMuonVetoPt<muonVetoPtThreshold) & ~(events.cscRechitClusterMuonVetoLooseId) &(events.cscRechitClusterJetVetoPt<30))

    if not isData:
        cluster_cuts = (cluster_cuts) & (events.cscRechitCluster_match_gLLP)
    else:
        cluster_cuts = (cluster_cuts) & (abs(events.cscRechitClusterMet_dPhi)>2)
    # if clusterSizeCut:
    #     cluster_cuts = (cluster_cuts) & (events.cscRechitClusterSize>=160)

    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, 
                local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
    
    DNN_Scores = events.cscRechitClusterDNN_bkgMC_plusBeamHalo[cluster_cuts]
    DNN_Scores = ak.flatten(DNN_Scores[ak.num(DNN_Scores)>0]).compute()

    cluster_size = events.cscRechitClusterSize[cluster_cuts]
    cluster_size = ak.flatten(cluster_size[ak.num(cluster_size)>0]).compute()

    cluster_z = events.cscRechitClusterZ[cluster_cuts]
    cluster_z = ak.flatten(cluster_z[ak.num(cluster_z)>0]).compute()
    
    client.close()

    return DNN_Scores, cluster_size, cluster_z