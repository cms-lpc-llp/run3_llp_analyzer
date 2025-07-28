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
#from dask_awkward.repartition import repartition_like

sys.path.append('.')
import Processing_Helpers

cfg_file_path = os.environ["CMSSW_BASE"] + "/src/run3_llp_analyzer/python/HNL_Plotting_HelperFunctions/cuts_config/"

signal_xSec = 4.255
processed_lumi = 109080
generated_signalEvents = 500000

signal_normalization_factor = processed_lumi*signal_xSec/generated_signalEvents


collection_masks = {
    "tau": [],
    "cscCluster": [],
    "event": []
}

def computeEfficiency(cumulative_events, Run, invert=False, **args):
    #print(args)
    #print(args["collection"])
    #print((((cumulative_events[args['branch']]>=args['minVal']) & (cumulative_events[args['branch']]<=args['maxVal'])).compute()))
    #print((args[args['collection']].compute()))
    #print(ak.type(cumulative_events[args['branch']]>=args['minVal']))
    #print(Run, args["collection"])
    
    if Run==2 and args["collection"]=="cscCluster":
        branch = args['branch']
        branch = branch.replace("cscRechitCluster", "cscRechitCluster3")
        single_event_mask = ((cumulative_events[branch]>=args['minVal']) & ((cumulative_events[branch]<=args['maxVal'])))
    else:
        single_event_mask = ((cumulative_events[args['branch']]>=args['minVal']) & ((cumulative_events[args['branch']]<=args['maxVal'])))
    if invert:
        single_event_mask = ~(single_event_mask)

    existing_mask = (args[args['collection']])

    #print(ak.all(ak.num(single_event_mask, axis=1)==ak.num(existing_mask)).compute())
    #single_event_mask = (single_event_mask) & (ak.ones_like(existing_mask, dtype=bool))
    #print(ak.type(single_event_mask))
    #print(ak.type(existing_mask))
    mask = (single_event_mask) & (existing_mask)
    #print("MASK", mask.compute())
    #print("updated mask type: ", ak.type(mask))
    if args['collection']!="event":
        mask_event = ak.any(mask,axis=1) # make it event mask if tau or cluster
    else:
        mask_event = mask
    pass_selection = cumulative_events[mask_event]
    #print("About to return mask.type: ", ak.type(mask))
    return pass_selection, len(pass_selection.evtNum.compute())/len(cumulative_events.evtNum.compute()), args['collection'], mask


def makeCutflow(events, cfg_file, isMC=False, noGenCuts=True, Run=3):
    '''
    code for producing cutflow for HNL+Tau analyzer
    takes dask-awkward array as input
    '''

    #apply basic noise filters preselection

    if isMC:
        cuts_efficiencies = pd.DataFrame(columns = ["Cut Name", "Number of Events", "Cut Efficiency", "Cumulative Efficiency"])
    else:
        cuts_efficiencies = pd.DataFrame(columns = ["Cut Name", "Number of Events", "Cut Efficiency", "Cumulative Efficiency"])

    if Run==3:
        events = events[(events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)]
    #else:
        #events = events[(events.Flag_all)] #& (events.Flag_ecalBadCalibFilter)]
    #print(len(events.evtNum.compute()))

    if isMC and noGenCuts==False:
        total_events_num = len(events.evtNum.compute())
        if Run==3:
            accepted_events = events[ak.any(events.gLLP_csc, axis=1)]
        else:
            #print(ak.values_astype(events.gLLP_csc, bool).compute())
            #accepted_events = events[ak.any(ak.values_astype(events.gLLP_csc, bool), axis=1)]
            #print(events.gLLP_csc[:,0].compute())
            accepted_events = events[events.gLLP_csc[:,0]==1]
            #print(accepted_events.evtNum.compute())
        accepted_events_num = len(accepted_events.evtNum.compute())
        
    else:
        accepted_events = events
        accepted_events_num = len(events.evtNum.compute())

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

    

    #print("checking dims", ak.all(ak.num(collection_masks['cscCluster'], axis=1)==ak.num(cumulative_events.cscRechitClusterSize, axis=1)).compute())
            
    # print(ak.type(cumulative_events.cscRechitClusterSize.compute()))
    # print(ak.type(collection_masks["cscCluster"].compute()))

    #load in selections
    if isMC:
        with open(cfg_file_path+"/"+cfg_file, 'r') as f:
            cuts_dict = yaml.safe_load(f)
        normalization_factor = signal_normalization_factor
    else:
        with open(cfg_file_path+"/"+cfg_file, 'r') as f:
            cuts_dict = yaml.safe_load(f)
        normalization_factor = 1

    print(normalization_factor)
    for cut, cut_info in cuts_dict.items():
        #print(cut)
        invert= "InnerRing" in cut
        # if cut=="pass_trigger" and Run==3:
        #     forward_mask = (cumulative_events.cscRechitClusterNRechitChamberMinus11==0) & (cumulative_events.cscRechitClusterNRechitChamberMinus12==0) & (cumulative_events.cscRechitClusterNRechitChamberPlus11==0) & (cumulative_events.cscRechitClusterNRechitChamberPlus12==0)
        #     forward_mask_total = ak.any((collection_masks['cscCluster']) & (forward_mask), axis=1)
        #     print("forward veto eff: ", len(cumulative_events.evtNum[forward_mask_total].compute())/len((cumulative_events.evtNum).compute()))
        if Run==2 and (cut_info["collection"]=="tau" or "tau" in cut or "trigger" in cut):
            continue
        cumulative_events, eff, collection, mask = computeEfficiency(cumulative_events, Run, invert, **cut_info, **collection_masks)
        #print("num events: ", len(cumulative_events.evtNum.compute()))
        #print("checking dims", ak.all(ak.num(collection_masks['cscCluster'][ak.any(mask, axis=1)], axis=1)==ak.num(cumulative_events.cscRechitClusterSize, axis=1)).compute())
        #print("checking dims", ak.all(ak.num(mask[ak.any(mask, axis=1)], axis=1)==ak.num(cumulative_events.cscRechitClusterSize, axis=1)).compute())
        #print("mask type: ", ak.type(mask))
        for key in list(collection_masks.keys()):
            #print(key, collection)
            #print(mask.compute(), collection_masks[key].compute())
            if key!=collection and collection!="event":
                collection_masks[key] = (collection_masks[key][ak.any(mask, axis=1)])
                #print("new mask len: ", len(collection_masks[key].compute()))
            elif collection!="event":
                existing_mod = collection_masks[key][ak.any(mask, axis=1)]
                new_mod = mask[ak.any(mask, axis=1)]
                collection_masks[key] = (existing_mod) & (new_mod)
                #print("checking dims", ak.all(ak.num(collection_masks[key], axis=1)==ak.num(cumulative_events.cscRechitClusterSize, axis=1)).compute())
                #print("new mask: ", collection_masks[key].compute())
            else:
                collection_masks[key] = collection_masks[key][mask]

        current_eff = eff*current_eff
        cuts_efficiencies.loc[len(cuts_efficiencies)] = [cut_info['name'], len(cumulative_events.evtNum.compute())*normalization_factor, f'{eff:.3f}', "%.3f"%current_eff]
        
    #topology cuts
    print("computing delta Eta cartesian")
    tau_cscCluster_eta = ak.cartesian([cumulative_events.tauEta[collection_masks['tau']], cumulative_events.cscRechitClusterEta[collection_masks['cscCluster']]], nested=True)
    eta_diffs = tau_cscCluster_eta["0"]-tau_cscCluster_eta["1"]
    #flatten properly, due to bug in this version of awkwaard
    eta_diffs = ak.Array([ak.concatenate([nested for nested in event], axis=0) for event in eta_diffs.compute()])
    single_mask = (eta_diffs>-2) & (eta_diffs<2)
    #single_mask = dak.from_awkward(single_mask, npartitions=1)
    
    single_mask_event = ak.any(single_mask, axis=1)
    
    existing_mask = collection_masks['event'].compute()
    
    mask = (single_mask_event) & (existing_mask)
    mask_event = mask
    
    eff = ak.count_nonzero(mask_event)/len(cumulative_events.evtNum.compute())
    #mask_dask = dak.from_awkward(mask_event, npartitions=cumulative_events.npartitions)
    evtNums = cumulative_events.evtNum.compute()[mask_event]
    tauPhi = cumulative_events.tauPhi.compute()[mask_event]
    clusterPhi = cumulative_events.cscRechitClusterPhi.compute()[mask_event]
    #print("Events len: ", len(cumulative_events.evtNum.compute()))
    current_eff = current_eff*eff
    collection="event"
    for key in list(collection_masks.keys()):
        print(key, collection)
        collection_masks[key] = collection_masks[key].compute()
        collection_masks[key] = collection_masks[key][mask_event]

    
    cuts_efficiencies.loc[len(cuts_efficiencies)] = ['|dEta(cluster, tau)<2|', len(evtNums)*normalization_factor, f'{eff:.3f}', "%.3f"%current_eff]
    

    #topology cuts
    print("computing delta Phi cartesian")
    tau_cscCluster_phi = ak.cartesian([tauPhi[collection_masks['tau']], clusterPhi[collection_masks['cscCluster']]], nested=True)
    phi_diffs = tau_cscCluster_phi["0"]-tau_cscCluster_phi["1"]
    phi_diffs = Processing_Helpers.deltaPhiAk(phi_diffs)
    #flatten properly, due to bug in this version of awkwaard
    phi_diffs = ak.Array([ak.concatenate([nested for nested in event], axis=0) for event in phi_diffs])
    single_mask = (phi_diffs>2) | (phi_diffs<-2)
    #single_mask = dak.from_awkward(single_mask, npartitions=1)
    
    single_mask_event = ak.any(single_mask, axis=1)
    
    existing_mask = collection_masks['event']
    
    mask = (single_mask_event) & (existing_mask)
    mask_event = mask
    
    eff = ak.count_nonzero(mask_event)/len(evtNums)
    evtNums = evtNums[mask_event]
    #mask_dask = dak.from_awkward(mask_event, npartitions=cumulative_events.npartitions)
    #cumulative_events = cumulative_events[mask_event]
    current_eff = current_eff*eff
    collection="event"
    for key in list(collection_masks.keys()):
        print(key, collection)
        #collection_masks[key] = collection_masks[key].compute()
        collection_masks[key] = collection_masks[key][mask_event]

    
    cuts_efficiencies.loc[len(cuts_efficiencies)] = ['|dPhi(cluster, tau)>2|', len(evtNums)*normalization_factor, f'{eff:.3f}', "%.3f"%current_eff]


    return cuts_efficiencies