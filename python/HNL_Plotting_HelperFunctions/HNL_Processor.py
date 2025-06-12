#!/usr/bin/env python3
import coffea
import awkward as ak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot

import hist

#helper function for processor class
#processor structure inspired by Martin Kwok https://github.com/kakwok/LLP_coffea/blob/main/HNLprocessor/HNLproc_4.py

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

def buildClusterrMask(event_mask, cluster_mask, tau_mask):
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


class HNL_Processor(processor.ProcessorABC):
    '''
    Code to proccess events dask-awkward arrays from output of HNL analyzers. Event info is factored into tau, cscRechitCluster, dtRechitCluster, and general event variables (noise filters, etc)
    '''
    
    def __init__(self,**options):
        defaultOptions = {'campaign':'MC_Summer24'}
        options = {**defaultOptions, **options}
        self.campaign = options["campaign"]
        self._accumulator = {"tauPt": hist.Hist(hist.axis.Regular(50,0,100, name="tauPt", label="Tau pT [GeV]", underflow=False, overflow=False))}

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
        event_mask = ak.zip({"nTau_mask":nTau})
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

        #generate analysis objects
        taus = self.buildRecoTaus(events)
        cscClusters = self.buildCscRechitClusters(events)
        
        #apply cuts to analysis objects and do event-level selections
        cscCluster_mask = self.cscClusterSelections(cscClusters)
        tau_mask = self.tauSelections(taus)
        event_mask = self.eventSelections(events)
        total_mask = buildTauMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask)
        print(ak.flatten((events.tauPt[total_mask]).compute()))
        output["tauPt"].fill(tauPt=ak.flatten(events.tauPt[total_mask].compute()))
        return output



    def postprocess(self, accumulator):
        return accumulator