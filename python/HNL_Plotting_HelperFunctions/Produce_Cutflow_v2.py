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
#from dask_awkward.repartition import repartition_like

sys.path.append('.')
import Processing_Helpers

cfg_file_path = os.environ["CMSSW_BASE"] + "/src/run3_llp_analyzer/python/HNL_Plotting_HelperFunctions/cuts_config/"

signal_xSec = 1
processed_lumi = 109080



collection_masks = {
    "tau": [],
    "cscCluster": [],
    "event": []
}

def computeEfficiency(cumulative_events, Run, **args):
    #print("invert id", args['invert'])
    invert = ('yes' in args['invert'] or 'Yes' in args['invert'])
    if Run==2 and args["collection"]=="cscCluster":
        branch = args['branch']
        branch = branch.replace("cscRechitCluster", "cscRechitCluster3")
        single_event_mask = ((cumulative_events[branch]>=args['minVal']) & ((cumulative_events[branch]<=args['maxVal'])))
    else:
        single_event_mask = ((cumulative_events[args['branch']]>=args['minVal']) & ((cumulative_events[args['branch']]<=args['maxVal'])))
    if invert:
        single_event_mask = ~(single_event_mask)

    existing_mask = (args[args['collection']])

    
    mask = (single_event_mask) & (existing_mask)
    
    if args['collection']!="event":
        mask_event = ak.any(mask,axis=1) # make it event mask if tau or cluster
    else:
        mask_event = mask
    pass_selection = cumulative_events[mask_event]
    
    return pass_selection, ak.sum(pass_selection.weights)/ak.sum(cumulative_events.weights), args['collection'], mask


def makeCutflow(events, cfg_file, isMC=False, noGenCuts=True, Run=3, sample_ctau=1000, reweight_ctau = 1000):
    '''
    code for producing cutflow for HNL+Tau analyzer
    takes dask-awkward array as input
    '''

    generated_signalEvents = ak.count(events.evtNum)
    print(generated_signalEvents.compute())
    signal_normalization_factor = processed_lumi*signal_xSec/generated_signalEvents
    print(signal_normalization_factor.compute())
    #print(len(events.evtNum.compute()))

    if isMC:
        weights = sample_ctau/reweight_ctau*np.exp(10*events["gLLP_ctau"]*(1/sample_ctau-1/reweight_ctau))
        #print(ak.sum(weights).compute())
        weights = weights/ak.sum(weights)*ak.count(events.evtNum)
        #make L1 weight event level - can't add to generic weight because branch is not filled if there is no cluster
        L1_weights_eventLevel = 1 - ak.prod(1 - events.cscRechitClusterHMTEfficiency, axis=1)
        #print(L1_weights_eventLevel.compute())
        
    else:
        weights = ak.ones_like(events["runNum"])
        L1_weights_eventLevel = ak.ones_like(events["runNum"])

    events = ak.with_field(events, L1_weights_eventLevel, "L1_weights_eventLevel")
    #print(weights.compute())
    #print(ak.sum(weights.compute()))
    events = ak.with_field(events, weights, "weights")
    
    #apply basic noise filters preselection

    cuts_efficiencies = pd.DataFrame(columns = ["Cut Name", "Number of Weighted Events", "Number of Unweighted Events", "Cut Efficiency", "Cumulative Efficiency"])
    #cuts_efficiencies = dd.from_pandas(cuts_efficiencies, npartitions=1)

    if Run==3:
        events = events[(events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)]
    
    if isMC and noGenCuts==False:
        if Run==3:
            accepted_events = events[ak.any(events.gLLP_csc, axis=1)]
        else:
            accepted_events = events[events.gLLP_csc[:,0]==1]
        
    else:
        accepted_events = events

    cumulative_events = accepted_events
    current_eff = 1

    #masks to keep track of WHICH tau or cluster is passing the cuts. Initial preselection is required to not have empty entries in the array
    if Run==2:
        collection_masks = {
        #"tau": cumulative_events.tauPt>0,
        "cscCluster":cumulative_events.cscRechitCluster3Size>0,
        "event": cumulative_events.evtNum>0
    }
    else:
        collection_masks = {
            "tau": cumulative_events.tauPt>0,
            "cscCluster":cumulative_events.cscRechitClusterSize>0,
            "event": cumulative_events.evtNum>0
        }

    #load in selections
    if isMC:
        with open(cfg_file_path+"/"+cfg_file, 'r') as f:
            cuts_dict = yaml.safe_load(f)
        normalization_factor = signal_normalization_factor
    else:
        with open(cfg_file_path+"/"+cfg_file, 'r') as f:
            cuts_dict = yaml.safe_load(f)
        normalization_factor = 1

    names, event_counts, event_counts_unweighted, effs, cumulative_effs = [],[],[],[], []
    for cut, cut_info in cuts_dict.items():
        
        #invert=False
        if Run==2 and (cut_info["collection"]=="tau" or "tau" in cut or "trigger" in cut):
            continue
        #print(cut_info)
        cumulative_events, eff, collection, mask = computeEfficiency(cumulative_events, Run, **cut_info, **collection_masks)
        
        for key in list(collection_masks.keys()):
            if key!=collection and collection!="event":
                collection_masks[key] = (collection_masks[key][ak.any(mask, axis=1)])
            elif collection!="event":
                existing_mod = collection_masks[key][ak.any(mask, axis=1)]
                new_mod = mask[ak.any(mask, axis=1)]
                collection_masks[key] = (existing_mod) & (new_mod)
            else:
                collection_masks[key] = collection_masks[key][mask]

        current_eff = eff*current_eff
        
        names.append(cut_info['name'])
        event_counts.append(ak.sum(cumulative_events.weights)*normalization_factor)
        event_counts_unweighted.append(ak.count(cumulative_events.weights))
        effs.append(eff)
        cumulative_effs.append(current_eff)
    
    client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
    
    print("about to compute")
    cuts_efficiencies = pd.DataFrame({"Cut Name":names, 
            "Number of Weighted Events":dask.compute(*event_counts), 
            "Number of Unweighted Events":dask.compute(*event_counts_unweighted), 
            "Cut Efficiency":dask.compute(*effs), 
            "Cumulative Efficiency":dask.compute(*cumulative_effs)
            })
    
    client.close()

    return cuts_efficiencies