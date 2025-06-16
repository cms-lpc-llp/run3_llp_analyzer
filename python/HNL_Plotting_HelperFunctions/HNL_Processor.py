#!/usr/bin/env python3

import os

import coffea
import awkward as ak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot

import hist
import yaml

#helper functions for processor class
#processor structure inspired by Martin Kwok https://github.com/kakwok/LLP_coffea/blob/main/HNLprocessor/HNLproc_4.py

def produce_hist_dict(cfg_file: str)->dict:
    '''
    return a dictionary with hist.Hist objects specfieid in cfg_file
    '''
    with open(cfg_file, 'r') as f:
        plot_cfgs = yaml.safe_load(f)
    hist_dict = {}
    for name, plot_info in plot_cfgs.items():
        hist_dict[name] = hist.Hist(hist.axis.Regular(plot_info['nbins'],plot_info['xmin'],plot_info['xmax'], 
            name="plot", label=plot_info['x_label'], underflow=plot_info['underflow'], 
            overflow=plot_info['overflow'],
            ),metadata={"title":plot_info["title"], "y_label":plot_info["y_label"],
                    "branch":plot_info["MuonSystem_Branch_Expression"]})

    return hist_dict


def combineMasks(masks):
    '''
    Take input of an awkward array where each field corresponds to a boolean array.
    It must be the same size for all fields, and the output is a single array with these dimensions
    '''
    mask = ak.ones_like(masks[masks.fields[0]])
    for field in masks.fields:
        mask = (mask) & (masks[field])
    return mask


def combineMasksFlatten(masks):
    '''
    Same as combine masks, but we take the AND along the object-level axis so the returned array consists of one bool per event
    '''
    flattened_mask = ak.all(ak.ones_like(masks[masks.fields[0]]),axis=1)
    for field in masks.fields:
        flattened_mask = (flattened_mask) & (ak.all(masks[field],axis=1))
    return flattened_mask

def buildEventMask(event_mask, cluster_mask, tau_mask):
    '''Build mask to apply to event-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasksFlatten(cluster_mask)
    total_tau_mask = combineMasksFlatten(tau_mask)
    total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)

    return total_mask

def buildClusterMask(event_mask, cluster_mask, tau_mask):
    '''Build mask to apply to event-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasks(cluster_mask)
    total_tau_mask = combineMasksFlatten(tau_mask)
    total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)

    return total_mask


def buildTauMask(event_mask, cluster_mask, tau_mask):
    '''Build mask to apply to event-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasksFlatten(cluster_mask)
    total_tau_mask = combineMasks(tau_mask)
    total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)

    return total_mask





######## Define processor class #############
#############################################

class HNL_Processor(processor.ProcessorABC):
    '''
    Code to proccess events dask-awkward arrays from output of HNL analyzers. Event info is factored into tau, cscRechitCluster, dtRechitCluster, and general event variables (noise filters, etc)
    '''
    path_to_configs = os.environ["CMSSW_BASE"]+"src/run3_llp_analyzer/python/HNL_Plotting_HelperFunctions/hist_configs/"

    
    def __init__(self,**options):
        defaultOptions = {'campaign':'MC_Summer24',
                          'tau_hists_config':'tau_hists.yaml',
                          'cscCluster_hists_config':'cscCluster_hists.yaml',
                          'eventLevel_hists_config':'eventLevel_hists.yaml'}
        options = {**defaultOptions, **options}
        
        self.campaign = options["campaign"]
        self.tau_hist_config = self.path_to_configs+options["tau_hists_config"]
        self.cscCluster_hist_config = self.path_to_configs+options["cscCluster_hists_config"]
        self.eventLevel_hists_config = self.path_to_configs+options["eventLevel_hists_config"]
        
        self.tau_hists_dict = produce_hist_dict(self.tau_hist_config)
        self.cscCluster_hists_dict = produce_hist_dict(self.cscCluster_hist_config)
        self.eventLevel_hists_dict = produce_hist_dict(self.eventLevel_hists_config)
        
        
        self._accumulator = {}

    @property
    def accumulator(self):
        return self._accumulator

    def buildCscRechitClusters(self,events):
        '''
        Build CSC Rechit Clusters from MuonSystem Output
        '''
        cscBranches = [branch for branch in events.fields if "tau" in branch and "ctau" not in branch]
        cscClusters = ak.zip({cscBranch:events[cscBranch] for cscBranch in cscBranches})
        return cscClusters

    def buildRecoTaus(self, events):
        '''
        Build reco taus from MuonSystem Output
        '''
        tauBranches = [branch for branch in events.fields if "tau" in branch and "ctau" not in branch]
        taus = ak.zip({tauBranch:events[tauBranch] for tauBranch in tauBranches})
        return taus

    
    def eventSelections(self, events):
        initial_mask = ak.ones_like(events, dtype=bool)
        nTau = events.nTaus==1
        noiseFilters = (events.Flag_all) & (events.jetVeto)
        event_mask = ak.zip({"nTau_mask":nTau, "noiseFilters":noiseFilters})
        return event_mask
    
    def cscClusterSelections(self, cscClusters):
        initial_mask = ak.ones_like(cscClusters, dtype=bool)
        return initial_mask

    def tauSelections(self, taus):
        initial_mask = ak.ones_like(taus, dtype=bool)
        tau_eta = abs(taus.tauEta)<2.3
        tau_mask = ak.zip({"eta_cut":tau_eta})
        return tau_mask



    def process(self, events):
        
        output = self.accumulator.copy()
        print(output)
        #generate analysis objects
        taus = self.buildRecoTaus(events)
        cscClusters = self.buildCscRechitClusters(events)
        
        #apply cuts to analysis objects and do event-level selections
        cscCluster_mask = self.cscClusterSelections(cscClusters)
        tau_mask = self.tauSelections(taus)
        event_mask = self.eventSelections(events)
        total_mask_tau = buildTauMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask)
        total_mask_cscCluster = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask)
        total_mask_event = buildEventMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask)
        
        for plot in self.tau_hists_dict.keys():
            tauHist = self.tau_hists_dict[plot]
            tauHist.fill(plot=ak.flatten(events[tauHist.metadata["branch"]][total_mask_tau].compute()))
            output[plot] = tauHist
        # for plot in self.cscCluster_hists_dict.keys():
        #     cscClusterHist = self.cscCluster_hists_dict[plot]
        #     cscClusterHist.fill(plot=ak.flatten(events[cscClusterHist.metadata["branch"]][total_mask_cscCluster].compute()))
        #     output[plot] = cscClusterHist
        for plot in self.eventLevel_hists_dict.keys():
            eventLevelHist = self.eventLevel_hists_dict[plot]
            eventLevelHist.fill(plot=events[eventLevelHist.metadata["branch"]][total_mask_event].compute())
            output[plot] = eventLevelHist
        return output



    def postprocess(self, accumulator):
        return accumulator