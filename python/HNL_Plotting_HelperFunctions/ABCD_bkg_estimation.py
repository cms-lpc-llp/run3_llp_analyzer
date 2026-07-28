#!/usr/bin/env python3

import sys
import os

import coffea
import awkward as ak
from coffea import processor

#from coffea.nanoevents.methods import vector
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


def run_ABCD_from_cutflow(events, cluster_mask, sizeCut, dPhiCut, normalization_factor=1, blind=True, flavor="Tau", isMC=False):
    '''
    Code to return dataframe with event in each bin and expected number of events in signal bin
    '''
    effective_size_cut = sizeCut
    if isMC:
        effective_size_cut = modeling_cut_lookup.remap_abcd_size_cut_for_mc(sizeCut, context="ABCD")
    print(f"[ABCD] Using sizeCut={effective_size_cut}, dPhiCut={dPhiCut}, isMC={isMC}")

    # Ensure each event contributes at most one cluster: pick max-size cluster per event
    sizes = events.cscRechitClusterSize
    dphi = abs(events[f'cscRechitClusterPrompt{flavor}DeltaPhi'])
    sizes_masked = ak.where(cluster_mask, sizes, -1)
    has_cluster = ak.any(cluster_mask, axis=1)
    idx = ak.argmax(sizes_masked, axis=1, keepdims=True)
    size_sel = ak.firsts(sizes[idx])
    dphi_sel = ak.firsts(dphi[idx])
    size_sel = ak.fill_none(size_sel, -1)
    dphi_sel = ak.fill_none(dphi_sel, 0)
    weights_sel = events.weights

    event_counts = [] #list to store dask objects (number of events in each bin) that need to be stored
    event_counts_unc = []
    if blind:
        bin_names = ["bin A (low NHits, low dPhi)", "bin B (low NHits, high dPhi)", "bin C (high NHits, low dPhi)", "bin D expected (high NHits, high dPhi)"]
    else:
        bin_names = ["bin A (low NHits, low dPhi)", "bin B (low NHits, high dPhi)", "bin C (high NHits, low dPhi)", "bin D (high NHits, high dPhi)"]
    
    #bin A
    mask_A = (has_cluster) & (size_sel<effective_size_cut) & (dphi_sel<dPhiCut)
    bin_A = ak.sum(weights_sel[mask_A])*normalization_factor
    #print(bin_A)
    event_counts.append(bin_A)
    bin_A_unc = (ak.sum(weights_sel[mask_A]**2)**0.5)*normalization_factor
    event_counts_unc.append(bin_A_unc)

    #bin B
    mask_B = (has_cluster) & (size_sel<effective_size_cut) & (dphi_sel>=dPhiCut)
    bin_B = ak.sum(weights_sel[mask_B])*normalization_factor
    print()
    event_counts.append(bin_B)
    bin_B_unc = (ak.sum(weights_sel[mask_B]**2)**0.5)*normalization_factor
    event_counts_unc.append(bin_B_unc)

    #bin C
    mask_C = (has_cluster) & (size_sel>=effective_size_cut) & (dphi_sel<dPhiCut)
    bin_C = ak.sum(weights_sel[mask_C])*normalization_factor
    event_counts.append(bin_C)
    bin_C_unc = (ak.sum(weights_sel[mask_C]**2)**0.5)*normalization_factor
    event_counts_unc.append(bin_C_unc)
    
    if blind:
    #bin D expected
        bin_D_exp = bin_B/bin_A*bin_C
        bin_D_exp_unc = bin_D_exp * ((bin_A_unc/bin_A)**2+(bin_B_unc/bin_B)**2+(bin_C_unc/bin_C)**2)**0.5
        event_counts.append(bin_D_exp)
        event_counts_unc.append(bin_D_exp_unc)
    else:
        mask_D = (has_cluster) & (size_sel>=effective_size_cut) & (dphi_sel>=dPhiCut)
        bin_D = ak.sum(weights_sel[mask_D])*normalization_factor
        bin_D_unc = (ak.sum(weights_sel[mask_D]**2)**0.5)*normalization_factor
        
        event_counts.append(bin_D)
        event_counts_unc.append(bin_D_unc)

    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
    
    event_counts_computed = dask.compute(*event_counts)
    events_counts_unc_computed = dask.compute(*event_counts_unc)

    #print(event_counts_computed)
    #print(events_counts_unc_computed)

    bin_counts_strs = [f"{event_counts_computed[i]} +- {events_counts_unc_computed[i]}" for i in range(len(event_counts))]

    client.close()
    return pd.DataFrame({"Bin":bin_names, "Counts":bin_counts_strs})

def run_ABCD_from_cutflow_dEta(events, cluster_mask, sizeCut, dEtaCut, normalization_factor=1, blind=True, flavor="Tau", isMC=False):
    '''
    Code to return dataframe with event in each bin and expected number of events in signal bin
    '''
    effective_size_cut = sizeCut
    if isMC:
        effective_size_cut = modeling_cut_lookup.remap_abcd_size_cut_for_mc(sizeCut, context="ABCD dEta")
    print(f"[ABCD] Using sizeCut={effective_size_cut}, dEtaCut={dEtaCut}, isMC={isMC}")

    # Ensure each event contributes at most one cluster: pick max-size cluster per event
    sizes = events.cscRechitClusterSize
    deta = abs(events[f'cscRechitClusterPrompt{flavor}DeltaEta'])
    sizes_masked = ak.where(cluster_mask, sizes, -1)
    has_cluster = ak.any(cluster_mask, axis=1)
    idx = ak.argmax(sizes_masked, axis=1, keepdims=True)
    size_sel = ak.firsts(sizes[idx])
    deta_sel = ak.firsts(deta[idx])
    size_sel = ak.fill_none(size_sel, -1)
    deta_sel = ak.fill_none(deta_sel, 0)
    weights_sel = events.weights

    event_counts = [] #list to store dask objects (number of events in each bin) that need to be stored
    event_counts_unc = []
    if blind:
        bin_names = ["bin A (low NHits, low dEta)", "bin B (low NHits, high dEta)", "bin C (high NHits, low dEta)", "bin D expected (high NHits, high dEta)"]
    else:
        bin_names = ["bin A (low NHits, low dEta)", "bin B (low NHits, high dEta)", "bin C (high NHits, low dEta)", "bin D (high NHits, high dEta)"]

    #bin A
    mask_A = (has_cluster) & (size_sel<effective_size_cut) & (deta_sel<dEtaCut)
    bin_A = ak.sum(weights_sel[mask_A])*normalization_factor
    #print(bin_A)
    event_counts.append(bin_A)
    bin_A_unc = (ak.sum(weights_sel[mask_A]**2)**0.5)*normalization_factor
    event_counts_unc.append(bin_A_unc)

    #bin B
    mask_B = (has_cluster) & (size_sel<effective_size_cut) & (deta_sel>=dEtaCut)
    bin_B = ak.sum(weights_sel[mask_B])*normalization_factor
    print()
    event_counts.append(bin_B)
    bin_B_unc = (ak.sum(weights_sel[mask_B]**2)**0.5)*normalization_factor
    event_counts_unc.append(bin_B_unc)

    #bin D
    mask_D = (has_cluster) & (size_sel>=effective_size_cut) & (deta_sel>=dEtaCut)
    bin_D = ak.sum(weights_sel[mask_D])*normalization_factor
    event_counts.append(bin_D)
    bin_D_unc = (ak.sum(weights_sel[mask_D]**2)**0.5)*normalization_factor
    event_counts_unc.append(bin_D_unc)

    if blind:
    #bin D expected
        bin_C_exp = bin_A/bin_B*bin_D
        bin_C_exp_unc = bin_C_exp * ((bin_A_unc/bin_A)**2+(bin_B_unc/bin_B)**2+(bin_D_unc/bin_D)**2)**0.5
        event_counts.append(bin_C_exp)
        event_counts_unc.append(bin_C_exp_unc)
    else:
        mask_C = (has_cluster) & (size_sel>=effective_size_cut) & (deta_sel<dEtaCut)
        bin_C = ak.sum(weights_sel[mask_C])*normalization_factor
        bin_C_unc = (ak.sum(weights_sel[mask_C]**2)**0.5)*normalization_factor
        
        event_counts.append(bin_C)
        event_counts_unc.append(bin_C_unc)

    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
    
    event_counts_computed = dask.compute(*event_counts)
    events_counts_unc_computed = dask.compute(*event_counts_unc)

    #print(event_counts_computed)
    #print(events_counts_unc_computed)

    bin_counts_strs = [f"{event_counts_computed[i]} +- {events_counts_unc_computed[i]}" for i in range(len(event_counts))]

    client.close()
    return pd.DataFrame({"Bin":bin_names, "Counts":bin_counts_strs})
