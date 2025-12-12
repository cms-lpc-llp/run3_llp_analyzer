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

sys.path.append('.')
import Processing_Helpers

#helper functions for processor class
#processor structure inspired by Martin Kwok https://github.com/kakwok/LLP_coffea/blob/main/HNLprocessor/HNLproc_4.py

collection_branch_mapping = {"tau": "nTaus", "cscRechitClusters": "nCscRechitClusters", "LLP": "nGLLP", "gTau": "nGenTau", "gVisTau":"nGenVisTau"}


def produce_hist_dict(cfg_file: str)->dict:
    '''
    return a dictionary with hist.Hist objects specfieid in cfg_file
    '''
    with open(cfg_file, 'r') as f:
        plot_cfgs = yaml.safe_load(f)
    
    hist_dict = {}
    for name, plot_info in plot_cfgs.items():
        if "mask_branch" not in list(plot_info.keys()):
            #initialize with dummy mask (should always pass)
            hist_dict[name] = hist.Hist(hist.axis.Regular(plot_info['nbins'],plot_info['xmin'],plot_info['xmax'], 
                name="plot", label=plot_info['x_label'], underflow=plot_info['underflow'], 
                overflow=plot_info['overflow'],
                ),metadata={"title":plot_info["title"], "y_label":plot_info["y_label"],
                        "branch":plot_info["MuonSystem_Branch_Expression"],
                        "mask_collection":"event", "mask_branch": "runNum", "mask_lowVal":0., "mask_highVal":np.inf})
        else:
            hist_dict[name] = hist.Hist(hist.axis.Regular(plot_info['nbins'],plot_info['xmin'],plot_info['xmax'], 
                name="plot", label=plot_info['x_label'], underflow=plot_info['underflow'], 
                overflow=plot_info['overflow'],
                ),metadata={"title":plot_info["title"], "y_label":plot_info["y_label"],
                        "branch":plot_info["MuonSystem_Branch_Expression"], 
                        "mask_collection": plot_info["mask_collection"], "mask_branch":plot_info["mask_branch"], "mask_lowVal":float(plot_info["mask_lowVal"]), "mask_highVal":float(plot_info["mask_highVal"])})
        
        
        if "Mask_to_Exclude" in list(plot_info.keys()):
            hist_dict[name].metadata['Mask_To_Exclude'] = plot_info['Mask_To_Exclude']
        else:
            hist_dict[name].metadata['Mask_To_Exclude'] = None

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
    Same as combine masks, but we take the OR along the object-level axis so the returned array consists of one bool per event
    Flatten at the very end, so that a single cluster/tau in the event has to pass all selections, and then flatten at the end so that it is an event level mask
    '''
    # flattened_mask = ak.any(ak.ones_like(masks[masks.fields[0]]),axis=1)
    # for field in masks.fields:
    #     flattened_mask = (flattened_mask) & (ak.any(masks[field],axis=1))
    mask = masks[masks.fields[0]]
    for field in masks.fields:
         mask = (mask) & (masks[field])
    flattened_mask = ak.any(mask, axis=1)
    return flattened_mask

def buildEventMask(event_mask, cluster_mask, tau_mask, gTauMask=None, gVisTauMask=None, gLLPMask=NanoEventsFactory, MC=True):
    '''Build mask to apply to event-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasksFlatten(cluster_mask)
    total_tau_mask = combineMasksFlatten(tau_mask)
    
    if MC:
        total_gTau_mask = combineMasksFlatten(gTauMask)
        total_gVisTau_mask = combineMasksFlatten(gVisTauMask)
        total_gLLP_mask = combineMasksFlatten(gLLPMask)
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)&(total_gTau_mask)&(total_gVisTau_mask)&(total_gLLP_mask)
    
    else:
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)
    
    return total_mask

def buildClusterMask(event_mask, cluster_mask, tau_mask, gTauMask=None, gVisTauMask=None, gLLPMask=None, MC=True):
    '''Build mask to apply to cluster-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasks(cluster_mask)
    total_tau_mask = combineMasksFlatten(tau_mask)

    if MC:
        total_gTau_mask = combineMasksFlatten(gTauMask)
        total_gVisTau_mask = combineMasksFlatten(gVisTauMask)
        total_gLLP_mask = combineMasksFlatten(gLLPMask)
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)&(total_gTau_mask)&(total_gVisTau_mask)&(total_gLLP_mask)

    else:
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)
    
    return total_mask


def buildTauMask(event_mask, cluster_mask, tau_mask, gTauMask=None, gVisTauMask=None, gLLPMask=None, MC=True):
    '''Build mask to apply to reco tau-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasksFlatten(cluster_mask)
    total_tau_mask = combineMasks(tau_mask)
    
    if MC:
        total_gTau_mask = combineMasksFlatten(gTauMask)
        total_gVisTau_mask = combineMasksFlatten(gVisTauMask)
        total_gLLP_mask = combineMasksFlatten(gLLPMask)
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)&(total_gTau_mask)&(total_gVisTau_mask)&(total_gLLP_mask)

    else:
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_tau_mask)
    
    return total_mask





######## Define processor class #############
#############################################

class HNL_Processor(processor.ProcessorABC):
    '''
    Code to proccess events dask-awkward arrays from output of HNL analyzers. Event info is factored into tau, cscRechitCluster, dtRechitCluster, and general event variables (noise filters, etc)
    '''
    path_to_configs = os.environ["CMSSW_BASE"] + "/src/run3_llp_analyzer/python/HNL_Plotting_HelperFunctions/hist_configs_tau/"

    
    def __init__(self,**options):
        defaultOptions = {'campaign':'MC_Summer24',
                          'tau_hists_config':'tau_hists.yaml',
                          'cscCluster_hists_config':'cscCluster_hists.yaml',
                          'eventLevel_hists_config':'eventLevel_hists.yaml', 
                          'gLLP_hists_config':'gLLP_hists.yaml',
                          'gTau_hists_config': 'genTau_hists.yaml',
                          'gVisTau_hists_config': 'genVisTau_hists.yaml',
                          'isMC':False,
                          'applyGenInfo':True,
                          'make_TauID_hists':False}
        options = {**defaultOptions, **options}
        
        self.campaign = options["campaign"]
        self.tau_hist_config = self.path_to_configs+options["tau_hists_config"]
        self.cscCluster_hist_config = self.path_to_configs+options["cscCluster_hists_config"]
        self.eventLevel_hists_config = self.path_to_configs+options["eventLevel_hists_config"]
        self.gLLP_hists_config  = self.path_to_configs+options["gLLP_hists_config"]
        self.gTau_hists_config = self.path_to_configs+options["gTau_hists_config"]
        self.gVisTau_hists_config = self.path_to_configs+options["gVisTau_hists_config"]
        self.isMC = options['isMC']
        self.applyGenInfo = options['applyGenInfo']
        self.make_TauID_hists = options['make_TauID_hists']
        self.tau_hists_dict = produce_hist_dict(self.tau_hist_config)
        self.cscCluster_hists_dict = produce_hist_dict(self.cscCluster_hist_config)
        self.eventLevel_hists_dict = produce_hist_dict(self.eventLevel_hists_config)
        self.gLLP_hists_dict = produce_hist_dict(self.gLLP_hists_config)
        self.gTau_hists_dict = produce_hist_dict(self.gTau_hists_config)
        self.gVisTau_hists_dict = produce_hist_dict(self.gVisTau_hists_config)
        
        self._accumulator = {}

    @property
    def accumulator(self):
        return self._accumulator

    def buildCscRechitClusters(self,events):
        '''
        Build CSC Rechit Clusters from MuonSystem Output
        '''
        cscBranches = [branch for branch in events.fields if "cscRechitCluster" in branch]
        cscClusters = ak.zip({cscBranch:events[cscBranch] for cscBranch in cscBranches})
        return cscClusters

    def buildRecoTaus(self, events):
        '''
        Build reco taus from MuonSystem Output
        '''
        tauBranches = [branch for branch in events.fields if "tau" in branch and "ctau" not in branch or 'deltaR_GenTauRecoTau' in branch]
        taus = ak.zip({tauBranch:events[tauBranch] for tauBranch in tauBranches})
        # tauID_WPs = ['VVVLoose', 'VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']
        # if tauID == None:
        #     initial_mask = ak.ones_like(taus.deltaR_GenTauRecoTau, dtype=bool)
        # elif tauID!=None and tauID not in tauID_WPs:
        #     raise ValueError("Specified Tau ID Working Point Not Recognized")
        # else:
        #     print("applying ID")
        #     initial_mask = taus['tauIs'+tauID]
        #     print(ak.count_nonzero(initial_mask).compute())
        return taus

    def buildGLLP(self, events):
        '''
        Build Gen LLP from MuonSystem Output
        '''
        gLLPBranches = [branch for branch in events.fields if "gLLP" in branch and ("cscRechitCluster" not in branch and "dtRechitCluster" not in branch)]
        gLLPs = ak.zip({gLLPBranch:events[gLLPBranch] for gLLPBranch in gLLPBranches})
        return gLLPs


    def buildGTaus(self, events):
        '''
        Build Gen Taus from MuonSystem Output
        '''
        gTauBranches = [branch for branch in events.fields if "gTau" in branch]
        gTaus = ak.zip({gTauBranch:events[gTauBranch] for gTauBranch in gTauBranches})
        return gTaus
    
    def buildGVisTaus(self, events):
        '''
        Build Gen Vis Taus from MuonSystem Output
        '''
        gVisTauBranches = [branch for branch in events.fields if "gVisTau" in branch and "DecayMode" not in branch]
        gVisTaus = ak.zip({gVisTauBranch:events[gVisTauBranch] for gVisTauBranch in gVisTauBranches})
        return gVisTaus

    
    def eventSelections(self, events, mask_to_exclude=None):
        noiseFilters = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
        triggerFilter = events.HLT_CscCluster100_PNetTauhPFJet10_Loose
        #triggerFilter = ak.ones_like(events.HLT_CscCluster100_PNetTauhPFJet10_Loose)
        nTausFilter = events.nTaus==1
        cscClusterFilter = events.nCscRechitClusters==1
        event_mask = ak.zip({"noiseFilters":noiseFilters, "triggerFilter":triggerFilter, "nTausFilter":nTausFilter, "cscClusterFilter":cscClusterFilter})
        #event_mask = ak.zip({"noiseFilters":noiseFilters, "triggerFilter":triggerFilter,"cscClusterFilter":cscClusterFilter})
        if self.isMC and self.applyGenInfo:
            event_mask = ak.with_field(event_mask, events.nGenVisTau>0, "nGenVisTau")
            #event_mask = ak.with_field(event_mask, ak.flatten(events.gTauMuDecay), "EMask")

        if not mask_to_exclude==None:
            event_mask = ak.without_field(event_mask, mask_to_exclude)

        return event_mask
    
    def cscClusterSelections(self, cscClusters, mask_to_exclude=None):
        initial_mask = ak.ones_like(cscClusters.cscRechitClusterSize, dtype=bool)
        cluster_mask = ak.zip({"initial_mask":initial_mask, "muonVeto": cscClusters.cscRechitClusterMuonVetoPt<30, "jetVeto": cscClusters.cscRechitClusterJetVetoPt<30, "cluster_size": cscClusters.cscRechitClusterSize>=100, "inTime": (cscClusters.cscRechitClusterTimeWeighted>-5) & (cscClusters.cscRechitClusterTimeWeighted<12.5)})
        #cluster_mask = ak.zip({"initial_mask":initial_mask})
        if self.isMC and self.applyGenInfo:
            cluster_mask = ak.with_field(cluster_mask, cscClusters.cscRechitCluster_match_gLLP, "match_gLLP")
            cluster_mask = ak.with_field(cluster_mask, (abs(cscClusters.cscRechitCluster_match_gLLP_decay_z)>400) & (abs(cscClusters.cscRechitCluster_match_gLLP_decay_z)<1100), "LLP_inCSCs")
        
        if not mask_to_exclude==None:
            cluster_mask = ak.without_field(cluster_mask, mask_to_exclude)
        
        return cluster_mask

    def tauSelections(self, taus, mask_to_exclude=None, invertTauId=False):
        
        initial_mask = ak.ones_like(taus.deltaR_GenTauRecoTau, dtype=bool)
        tau_ID_mask = taus.tauIsLoose==(not invertTauId) #change back when not looking at the cluster size!
        tau_pT_mask = taus.tauPt>18
        tau_mask = ak.zip({"tau_ID_mask":tau_ID_mask, "tau_pT_mask":tau_pT_mask, "tau_eta": abs(taus.tauEta)<2.5})
        if self.isMC and self.applyGenInfo:
            deltaR_cut= taus.deltaR_GenTauRecoTau<0.4
            tau_mask = ak.with_field(tau_mask, deltaR_cut, "deltaR_cut")
        
        if not mask_to_exclude==None:
            tau_mask = ak.without_field(tau_mask, mask_to_exclude)

        return tau_mask

    def gTauSelections(self, gTaus):
        return ak.zip({"eta_cut": abs(gTaus.gTauEta)<2.5})
    
    def gVisSelections(self, gVisTaus):
        return ak.zip({"eta_cut": abs(gVisTaus.gVisTauEta)<2.3}) #change back to 2.5 after ID reconstruction studies
    
    def gLLPSelections(self, gLLPs):
        return ak.zip({"no_cut": ak.ones_like(gLLPs.gLLP_pt, dtype=bool)})

    
    
    
    
    
    ###     PROCESSOR   ####
    ########################
    def process(self, events, hists_to_process: list=None, fillGenHists = True, tauID: str='', tau_cluster_topo_hists=True, branchNames_for_invertTauID = ["CSC_Cluster_Size"]):
        #tauID argument explicitly for ID reconstruction studies, NOT generic cutflow
        #events = events[events.HLT_CscCluster100_PNetTauhPFJet10_Loose]
        hist_list = []
        if hists_to_process==None and self.isMC: #this means process all hists, since no specific ones are passed to the processor
            hist_list = [*list(self.eventLevel_hists_dict.keys()), *list(self.cscCluster_hists_dict.keys()), *list(self.gLLP_hists_dict.keys()), *list(self.gTau_hists_dict.keys()), *list(self.gVisTau_hists_dict.keys()), *list(self.tau_hists_dict.keys())]
        elif hists_to_process==None:
            hist_list = [*list(self.eventLevel_hists_dict.keys()), *list(self.cscCluster_hists_dict.keys()), *list(self.tau_hists_dict.keys())]
        else:
            hist_list = hists_to_process
        if not self.make_TauID_hists:
            hist_list = [hist for hist in hist_list if "ID" not in hist]
        print(hist_list)
        output = self.accumulator.copy()
        
        #generate analysis objects
        taus = self.buildRecoTaus(events)
        cscClusters = self.buildCscRechitClusters(events)
        
        #apply cuts to analysis objects and do event-level selections
        cscCluster_mask = self.cscClusterSelections(cscClusters)
        tau_mask = self.tauSelections(taus)
        tau_mask_invertId = self.tauSelections(taus, invertTauId=True)
        event_mask = self.eventSelections(events)

        if self.isMC and self.applyGenInfo:
            gTaus = self.buildGTaus(events)
            gVisTaus = self.buildGVisTaus(events)
            gLLP = self.buildGLLP(events)
            gLLP_mask = self.gLLPSelections(gLLP)
            gTau_mask = self.gTauSelections(gTaus)
            gVisTau_mask = self.gVisSelections(gVisTaus)
            
            total_mask_tau = buildTauMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask, gLLPMask=gLLP_mask, gVisTauMask=gVisTau_mask, gTauMask=gTau_mask, MC=True)
            #change back cluster mask
            total_mask_cscCluster = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask, gLLPMask=gLLP_mask, gVisTauMask=gVisTau_mask, gTauMask=gTau_mask, MC=True)
            total_mask_cscCluster_invertTauId = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask_invertId, gLLPMask=gLLP_mask, gVisTauMask=gVisTau_mask, gTauMask=gTau_mask, MC=True)
            total_mask_event = buildEventMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask, gLLPMask=gLLP_mask, gVisTauMask=gVisTau_mask, gTauMask=gTau_mask, MC=True)
            total_mask_gLLP = combineMasks(gLLP_mask)
            total_mask_gTau = combineMasks(gTau_mask)
            total_mask_gVisTau = combineMasks(gVisTau_mask)
            
        
        else:
            #use dummy true mask for gen-level masks
            total_mask_gLLP = None
            total_mask_gTau = None
            total_mask_gVisTau = None
            total_mask_tau = buildTauMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask, MC=False)
            total_mask_cscCluster_invertTauId = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask_invertId, MC=False)
            total_mask_cscCluster = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask, MC=False)
            total_mask_event = buildEventMask(event_mask=event_mask, cluster_mask=cscCluster_mask, tau_mask=tau_mask, MC=False)

        #print(total_mask_event.compute())
        print("Generated Masks, Starting to Fill Event-Level Histograms")

        collection_mask_mapping = {"tau": self.tauSelections, "cscRechitCluster": self.cscClusterSelections, "event": self.eventSelections}
        collection_objects = {"tau": taus, "cscRechitCluster": cscClusters, "event": events}
        collection_totalMaskFunction_mapping = {"tau": buildTauMask, "cscCluster": buildClusterMask, "event": buildEventMask}
        collection_mask_dict = {"tau": total_mask_tau, "cscCluster": total_mask_cscCluster, "event": total_mask_event,
                                "gLLP": total_mask_gLLP, "gTau": total_mask_gTau, "gVisTau": total_mask_gVisTau}

        def computeModMask(excluded_mask_collection, excluded_mask):
            mod_mask = collection_mask_mapping[excluded_mask_collection](collection_objects[excluded_mask_collection], excluded_mask)
            mask_dict_tmp = copy.deepcopy(collection_mask_dict[excluded_mask_collection])
            mask_dict_tmp[excluded_mask_collection] = mask_dict_tmp
            event_mask = collection_totalMaskFunction_mapping[excluded_mask_collection](**mask_dict_tmp)
            return event_mask

        for plot in self.eventLevel_hists_dict.keys():
            #print(plot)
            if plot not in hist_list:continue
            eventLevelHist = self.eventLevel_hists_dict[plot]
            added_mask = (events[eventLevelHist.metadata["mask_branch"]]>=eventLevelHist.metadata["mask_lowVal"]) & (events[eventLevelHist.metadata["mask_branch"]]<=eventLevelHist.metadata["mask_highVal"])
            if eventLevelHist.metadata["mask_collection"]!="event":
                added_mask = (events[collection_branch_mapping[eventLevelHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
            if eventLevelHist.metadata['Mask_To_Exclude']!=None:
                event_mask = computeModMask(eventLevelHist.metadata['Mask_To_Exclude_Collection'], eventLevelHist.metadata['Mask_To_Exclude'])
            else:
                event_mask = total_mask_event
            hist_mask = (added_mask) & (event_mask)
            eventLevelHist.fill(plot=events[eventLevelHist.metadata["branch"]][hist_mask].compute())
            output[plot] = eventLevelHist

        if self.isMC and fillGenHists and self.applyGenInfo:
            print("filling gen-level hists")
            for plot in self.gLLP_hists_dict.keys():
                if plot not in hist_list:continue
                gLLPHist = self.gLLP_hists_dict[plot]
                added_mask = (events[gLLPHist.metadata["mask_branch"]]>=gLLPHist.metadata["mask_lowVal"]) & (events[gLLPHist.metadata["mask_branch"]]<=gLLPHist.metadata["mask_highVal"])
                if gLLPHist.metadata["mask_collection"]!="event":
                    added_mask = (events[collection_branch_mapping[gLLPHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
                hist_mask = (added_mask) & (total_mask_gLLP)
                gLLPHist.fill(plot=ak.flatten(events[gLLPHist.metadata["branch"]][hist_mask].compute()))
                output[plot] = gLLPHist
            
            for plot in self.gTau_hists_dict.keys():
                if plot not in hist_list:continue
                gTauHist = self.gTau_hists_dict[plot]
                added_mask = (events[gTauHist.metadata["mask_branch"]]>=gTauHist.metadata["mask_lowVal"]) & (events[gTauHist.metadata["mask_branch"]]<=gTauHist.metadata["mask_highVal"])
                if gTauHist.metadata["mask_collection"]!="event":
                    (events[collection_branch_mapping[gTauHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
                hist_mask = (added_mask) & (total_mask_gTau)
                gTauHist.fill(plot=ak.flatten(events[gTauHist.metadata["branch"]][hist_mask].compute()))
                output[plot] = gTauHist
            
            for plot in self.gVisTau_hists_dict.keys():
                if plot not in hist_list:continue
                gVisTauHist = self.gVisTau_hists_dict[plot]
                total_mask = ak.ones_like(events.evtNum,dtype=bool)
                added_mask = (events[gVisTauHist.metadata["mask_branch"]]>=gVisTauHist.metadata["mask_lowVal"]) & (events[gVisTauHist.metadata["mask_branch"]]<=gVisTauHist.metadata["mask_highVal"])
                if gVisTauHist.metadata["mask_collection"]!="event":
                    added_mask = (events[collection_branch_mapping[gVisTauHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
                    if gVisTauHist.metadata["mask_collection"]=="tau" and tauID!='': #FOR ID RECO STUDES
                        added_mask = (added_mask) & (combineMasksFlatten(ak.zip({"ID": events['tauIs'+tauID]})))
                hist_mask = (added_mask) & (total_mask_gVisTau)
                gVisTauHist.fill(plot=ak.flatten(events[gVisTauHist.metadata["branch"]][hist_mask].compute()))
                output[plot] = gVisTauHist
        
        print("filling reco taus hists")
        for plot in self.tau_hists_dict.keys():
            if plot not in hist_list:continue
            tauHist = self.tau_hists_dict[plot]
            added_mask = (events[tauHist.metadata["mask_branch"]]>=tauHist.metadata["mask_lowVal"]) & (events[tauHist.metadata["mask_branch"]]<=tauHist.metadata["mask_highVal"])
            if tauHist.metadata["mask_collection"]!="event":
                added_mask = (events[collection_branch_mapping[tauHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
            hist_mask = (added_mask) & (total_mask_tau)
            tauHist.fill(plot=ak.flatten(events[tauHist.metadata["branch"]][hist_mask].compute()))
            output[plot] = tauHist
        
        for plot in self.cscCluster_hists_dict.keys():
            if plot not in hist_list:continue
            print(plot)
            cscClusterHist = self.cscCluster_hists_dict[plot]
            added_mask = (events[cscClusterHist.metadata["mask_branch"]]>=cscClusterHist.metadata["mask_lowVal"]) & (events[cscClusterHist.metadata["mask_branch"]]<=cscClusterHist.metadata["mask_highVal"])
            if cscClusterHist.metadata["mask_collection"]!="event":
                added_mask = (events[collection_branch_mapping[cscClusterHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
            hist_mask = (added_mask) & (total_mask_cscCluster)
            if cscClusterHist.metadata['Mask_To_Exclude']!=None:
                print("excluding mask")
                cluster_mask = computeModMask(cscClusterHist.metadata['Mask_To_Exclude_Collection'], cscClusterHist.metadata['Mask_To_Exclude'])
            elif plot=="CSC_Cluster_Size" and not self.isMC:
                print("Inverting Tau LooseID in data for cluster size")
                cluster_mask = total_mask_cscCluster_invertTauId
            else:
                cluster_mask = total_mask_cscCluster
            hist_mask = (added_mask) & (cluster_mask)
            #hist_mask = (added_mask) & (combineMasks(cscCluster_mask)) & (combineMasks(event_mask)) # need to add back tau mask when I better understand where we're losing events
            cscClusterHist.fill(plot=ak.flatten(events[cscClusterHist.metadata["branch"]][hist_mask].compute()))
            output[plot] = cscClusterHist
        
    
        #filling histograms with delta phi/eta between tau and cscCluster - not obvious to me how to do it via config files
        
        
        if tau_cluster_topo_hists:
            print("filling topo hists")
            tau_cluster_dEta_hist = hist.Hist(hist.axis.Regular(30 ,-6,6, 
                name="plot", label="d$\eta$(tau, cluster)", underflow=False, 
                overflow=False,
                ),metadata={"title":"dEta(tau, cluster)", "y_label":"counts"})

            print("about to run cartesian for eta")
            tau_cscCluster_eta = ak.flatten(ak.flatten(ak.cartesian([events.tauEta[total_mask_tau], events.cscRechitClusterEta[total_mask_cscCluster]], nested=True)))
            tau_cluster_dEta_hist.fill(plot=(tau_cscCluster_eta["0"]-tau_cscCluster_eta["1"]).compute())
            output['tau_cluster_dEta'] = tau_cluster_dEta_hist

            print("about to run cartesian for phi")
            tau_cluster_dPhi_hist = hist.Hist(hist.axis.Regular(30 ,-6,6, 
                    name="plot", label="d$\phi$(tau, cluster)", underflow=False, 
                    overflow=False,
                    ),metadata={"title":"dPhi(tau, cluster)", "y_label":"counts"})



        
            tau_cscCluster_phi = ak.flatten(ak.flatten(ak.cartesian([events.tauPhi[total_mask_tau], events.cscRechitClusterPhi[total_mask_cscCluster]], nested=True)))
            tau_cscCluster_phi_diffs = Processing_Helpers.deltaPhi(tau_cscCluster_phi["0"]-tau_cscCluster_phi["1"])
            
            tau_cluster_dEta_hist.fill(plot=(tau_cscCluster_eta["0"]-tau_cscCluster_eta["1"]).compute())
            tau_cluster_dPhi_hist.fill(plot=tau_cscCluster_phi_diffs.compute())
            output['tau_cluster_dEta'] = tau_cluster_dEta_hist
            output['tau_cluster_dPhi'] = tau_cluster_dPhi_hist
        
        return output



    def postprocess(self, accumulator):
        return accumulator
