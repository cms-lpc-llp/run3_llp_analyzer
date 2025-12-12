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


def run_ABCD_from_cutflow(events, cluster_mask, sizeCut, dPhiCut, normalization_factor=1, blind=True, flavor="Tau"):
    '''
    Code to return dataframe with event in each bin and expected number of events in signal bin
    '''

    event_counts = [] #list to store dask objects (number of events in each bin) that need to be stored
    event_counts_unc = []
    if blind:
        bin_names = ["bin A (low NHits, low dPhi)", "bin B (low NHits, high dPhi)", "bin C (high NHits, low dPhi)", "bin D expected (high NHits, high dPhi)"]
    else:
        bin_names = ["bin A (low NHits, low dPhi)", "bin B (low NHits, high dPhi)", "bin C (high NHits, low dPhi)", "bin D (high NHits, high dPhi)"]
    
    #bin A
    bin_A = ak.sum(events.weights[ak.any((events.cscRechitClusterSize<sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaPhi'])<dPhiCut) & (cluster_mask), axis=1)])*normalization_factor
    #print(bin_A)
    event_counts.append(bin_A)
    bin_A_unc = bin_A**0.5
    event_counts_unc.append(bin_A_unc)

    #bin B
    bin_B = ak.sum(events.weights[ak.any((events.cscRechitClusterSize<sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaPhi'])>=dPhiCut) & (cluster_mask), axis=1)])*normalization_factor
    print()
    event_counts.append(bin_B)
    bin_B_unc = bin_B**0.5
    event_counts_unc.append(bin_B_unc)

    #bin C
    bin_C = ak.sum(events.weights[ak.any((events.cscRechitClusterSize>=sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaPhi'])<dPhiCut) & (cluster_mask), axis=1)])*normalization_factor
    event_counts.append(bin_C)
    bin_C_unc = bin_C**0.5
    event_counts_unc.append(bin_C_unc)
    
    if blind:
    #bin D expected
        bin_D_exp = bin_B/bin_A*bin_C
        bin_D_exp_unc = bin_D_exp * ((bin_A_unc/bin_A)**2+(bin_B_unc/bin_B)**2+(bin_C_unc/bin_C)**2)**0.5
        event_counts.append(bin_D_exp)
        event_counts_unc.append(bin_D_exp_unc)
    else:
        bin_D = ak.sum(events.weights[ak.any((events.cscRechitClusterSize>=sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaPhi'])>=dPhiCut) & (cluster_mask), axis=1)])*normalization_factor
        bin_D_unc = bin_D**0.5
        
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

def run_ABCD_from_cutflow_dEta(events, cluster_mask, sizeCut, dEtaCut, normalization_factor=1, blind=True, flavor="Tau"):
    '''
    Code to return dataframe with event in each bin and expected number of events in signal bin
    '''

    event_counts = [] #list to store dask objects (number of events in each bin) that need to be stored
    event_counts_unc = []
    if blind:
        bin_names = ["bin A (low NHits, low dEta)", "bin B (low NHits, high dEta)", "bin C (high NHits, low dEta)", "bin D expected (high NHits, high dEta)"]
    else:
        bin_names = ["bin A (low NHits, low dEta)", "bin B (low NHits, high dEta)", "bin C (high NHits, low dEta)", "bin D (high NHits, high dEta)"]

    #bin A
    bin_A = ak.sum(events.weights[ak.any((events.cscRechitClusterSize<sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaEta'])<dEtaCut) & (cluster_mask), axis=1)])*normalization_factor
    #print(bin_A)
    event_counts.append(bin_A)
    bin_A_unc = bin_A**0.5
    event_counts_unc.append(bin_A_unc)

    #bin B
    bin_B = ak.sum(events.weights[ak.any((events.cscRechitClusterSize<sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaEta'])>=dEtaCut) & (cluster_mask), axis=1)])*normalization_factor
    print()
    event_counts.append(bin_B)
    bin_B_unc = bin_B**0.5
    event_counts_unc.append(bin_B_unc)

    #bin D
    bin_D = ak.sum(events.weights[ak.any((events.cscRechitClusterSize>=sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaEta'])>=dEtaCut) & (cluster_mask), axis=1)])*normalization_factor
    event_counts.append(bin_D)
    bin_D_unc = bin_D**0.5
    event_counts_unc.append(bin_D_unc)

    if blind:
    #bin D expected
        bin_C_exp = bin_A/bin_B*bin_D
        bin_C_exp_unc = bin_C_exp * ((bin_A_unc/bin_A)**2+(bin_B_unc/bin_B)**2+(bin_D_unc/bin_D)**2)**0.5
        event_counts.append(bin_C_exp)
        event_counts_unc.append(bin_C_exp_unc)
    else:
        bin_C = ak.sum(events.weights[ak.any((events.cscRechitClusterSize>=sizeCut) & (abs(events[f'cscRechitClusterPrompt{flavor}DeltaEta'])<dEtaCut) & (cluster_mask), axis=1)])*normalization_factor
        bin_C_unc = bin_C**0.5
        
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

