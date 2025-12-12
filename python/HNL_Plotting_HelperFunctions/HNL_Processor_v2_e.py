#!/usr/bin/env python3

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
import sys
import copy
import dask
from dask.distributed import Client
import gc

sys.path.append('.')
import Processing_Helpers

#helper functions for processor class
#processor structure inspired by Martin Kwok https://github.com/kakwok/LLP_coffea/blob/main/HNLprocessor/HNLproc_4.py

collection_branch_mapping = {"electron": "nLeptons", "cscRechitClusters": "nCscRechitClusters", "LLP": "nGLLP", "gLep": "nGenLeptons"}
electron_pdgId = 11
electron_triggerPath = "HLT_CscCluster100_Ele5"


def produce_hist_dict(cfg_file: str, bin_multiplier: int=1)->dict:
    '''
    return a dictionary with hist.Hist objects specfieid in cfg_file
    '''
    #print(bin_multiplier)
    with open(cfg_file, 'r') as f:
        plot_cfgs = yaml.safe_load(f)
    
    hist_dict = {}
    for name, plot_info in plot_cfgs.items():
        if "mask_branch" not in list(plot_info.keys()):
            #initialize with dummy mask (should always pass)
            #print(type(plot_info['nbins']))
            hist_dict[name] = hist.Hist(hist.axis.Regular(int(plot_info['nbins'])*bin_multiplier,plot_info['xmin'],plot_info['xmax'], 
                name="plot", label=plot_info['x_label'], underflow=plot_info['underflow'], 
                overflow=plot_info['overflow'],
                ),metadata={"title":plot_info["title"], "x_label":plot_info["x_label"], "y_label":plot_info["y_label"],
                        "branch":plot_info["MuonSystem_Branch_Expression"],
                        "mask_collection":"event", "mask_branch": "runNum", "mask_lowVal":0., "mask_highVal":np.inf})
            #print(hist_dict[name].axes[0].size)
            #print(hist_dict[name].size)
        else:
            hist_dict[name] = hist.Hist(hist.axis.Regular(int(plot_info['nbins'])*bin_multiplier,plot_info['xmin'],plot_info['xmax'], 
                name="plot", label=plot_info['x_label'], underflow=plot_info['underflow'], 
                overflow=plot_info['overflow'],
                ),metadata={"title":plot_info["title"], "x_label":plot_info["x_label"], "y_label":plot_info["y_label"],
                        "branch":plot_info["MuonSystem_Branch_Expression"], 
                        "mask_collection": plot_info["mask_collection"], "mask_branch": plot_info["mask_branch"], "mask_lowVal":float(plot_info["mask_lowVal"]), "mask_highVal":float(plot_info["mask_highVal"])})
        
        
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
    Flatten at the very end, so that a single cluster/lepton in the event has to pass all selections, and then flatten at the end so that it is an event level mask
    '''
    # flattened_mask = ak.any(ak.ones_like(masks[masks.fields[0]]),axis=1)
    # for field in masks.fields:
    #     flattened_mask = (flattened_mask) & (ak.any(masks[field],axis=1))
    mask = masks[masks.fields[0]]
    for field in masks.fields:
         mask = (mask) & (masks[field])
    flattened_mask = ak.any(mask, axis=1)
    return flattened_mask

def buildEventMask(event_mask, cluster_mask, electron_mask, gLeptonMask=None, gLLPMask=None, MC=True):
    '''Build mask to apply to event-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasksFlatten(cluster_mask)
    total_electron_mask = combineMasksFlatten(electron_mask)
    
    if MC:
        total_gLepton_mask = combineMasksFlatten(gLeptonMask)
        total_gLLP_mask = combineMasksFlatten(gLLPMask)
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_electron_mask)&(total_gLepton_mask)&(total_gLLP_mask)
    
    else:
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_electron_mask)
    
    return total_mask

def buildClusterMask(event_mask, cluster_mask, electron_mask, gLeptonMask=None, gLLPMask=None, MC=True):
    '''Build mask to apply to cluster-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasks(cluster_mask)
    total_electron_mask = combineMasksFlatten(electron_mask)

    if MC:
        total_gLepton_mask = combineMasksFlatten(gLeptonMask)
        total_gLLP_mask = combineMasksFlatten(gLLPMask)
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_electron_mask)&(total_gLepton_mask)&(total_gLLP_mask)

    else:
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_electron_mask)
    
    return total_mask


def buildElectronMask(event_mask, cluster_mask, electron_mask, gLeptonMask=None, gLLPMask=None, MC=True):
    '''Build mask to apply to reco electron-level quantities'''
    
    total_event_mask = combineMasks(event_mask)
    total_cluster_mask = combineMasksFlatten(cluster_mask)
    total_electron_mask = combineMasks(electron_mask)
    
    if MC:
        total_gLepton_mask = combineMasksFlatten(gLeptonMask)
        total_gLLP_mask = combineMasksFlatten(gLLPMask)
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_electron_mask)&(total_gLepton_mask)&(total_gLLP_mask)

    else:
        total_mask = (total_event_mask)&(total_cluster_mask)&(total_electron_mask)
    
    return total_mask





######## Define processor class #############
#############################################

class HNL_Processor_v2_e(processor.ProcessorABC):
    '''
    Code to proccess events dask-awkward arrays from output of HNL analyzers. Event info is factored into lepton, cscRechitCluster, dtRechitCluster, and general event variables (noise filters, etc)
    '''
    path_to_configs = os.environ["CMSSW_BASE"] + "/src/run3_llp_analyzer/python/HNL_Plotting_HelperFunctions/hist_configs_electron/"

    
    def __init__(self,**options):
        defaultOptions = {'campaign':'MC_Summer24',
                          'electron_hists_config':'electron_hists.yaml',
                          'cscCluster_hists_config':'cscCluster_hists.yaml',
                          'eventLevel_hists_config':'eventLevel_hists.yaml', 
                          'gLLP_hists_config':'gLLP_hists.yaml',
                          'gLepton_hists_config': 'genElectron_hists.yaml',
                          'isMC':False,
                          'applyGenInfo':True,
                          'make_ElectronID_hists':False,
                          'bin_multiplier':1,
                          }
        options = {**defaultOptions, **options}
        
        self.campaign = options["campaign"]
        self.electron_hist_config = self.path_to_configs+options["electron_hists_config"]
        self.cscCluster_hist_config = self.path_to_configs+options["cscCluster_hists_config"]
        self.eventLevel_hists_config = self.path_to_configs+options["eventLevel_hists_config"]
        self.gLLP_hists_config  = self.path_to_configs+options["gLLP_hists_config"]
        self.gLepton_hists_config = self.path_to_configs+options["gLepton_hists_config"]
        self.isMC = options['isMC']
        self.applyGenInfo = options['applyGenInfo']
        self.make_ElectronID_hists = options['make_ElectronID_hists']
        self.electron_hists_dict = produce_hist_dict(self.electron_hist_config, bin_multiplier=options['bin_multiplier'])
        self.cscCluster_hists_dict = produce_hist_dict(self.cscCluster_hist_config, bin_multiplier=options['bin_multiplier'])
        self.eventLevel_hists_dict = produce_hist_dict(self.eventLevel_hists_config, bin_multiplier=options['bin_multiplier'])
        self.gLLP_hists_dict = produce_hist_dict(self.gLLP_hists_config, bin_multiplier=options['bin_multiplier'])
        self.gLepton_hists_dict = produce_hist_dict(self.gLepton_hists_config, bin_multiplier=options['bin_multiplier'])
        
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

    def buildRecoElectrons(self, events):
        '''
        Build reco electrons from MuonSystem Output
        '''
        electronBranches = [branch for branch in events.fields if "lep" in branch and "ctau" not in branch]
        #print(electronBranches)
        #print(ak.count_nonzero(abs(events.lepPdgId)==electron_pdgId).compute())
        electrons = ak.zip({electronBranch:events[electronBranch][abs(events.lepPdgId)==electron_pdgId] for electronBranch in electronBranches})

        # tauID_WPs = ['VVVLoose', 'VVLoose', 'VLoose', 'Loose', 'Medium', 'Tight', 'VTight', 'VVTight']
        # if tauID == None:
        #     initial_mask = ak.ones_like(taus.deltaR_GenTauRecoTau, dtype=bool)
        # elif tauID!=None and tauID not in tauID_WPs:
        #     raise ValueError("Specified Tau ID Working Point Not Recognized")
        # else:
        #     print("applying ID")
        #     initial_mask = taus['tauIs'+tauID]
        #     print(ak.count_nonzero(initial_mask).compute())
        return electrons

    def buildGLLP(self, events):
        '''
        Build Gen LLP from MuonSystem Output
        '''
        gLLPBranches = [branch for branch in events.fields if "gLLP" in branch and ("cscRechitCluster" not in branch and "dtRechitCluster" not in branch)]
        gLLPs = ak.zip({gLLPBranch:events[gLLPBranch] for gLLPBranch in gLLPBranches})
        return gLLPs


    def buildGLeptons(self, events):
        '''
        Build Gen Leptons from MuonSystem Output
        '''
        gLeptonBranches = [branch for branch in events.fields if "gLep" in branch]
        gLeptons = ak.zip({gLepBranch:events[gLepBranch] for gLepBranch in gLeptonBranches})
        return gLeptons
    
    def eventSelections(self, events, fromTau = False, mask_to_exclude=None):
        noiseFilters = (events.Flag_all) & (events.jetVeto) & (events.Flag_ecalBadCalibFilter)
        #noiseFilters = (events.Flag_all) & (events.Flag_ecalBadCalibFilter)
        triggerFilter = events[electron_triggerPath]
        #triggerFilter = ak.ones_like(events.HLT_CscCluster100_PNetTauhPFJet10_Loose)
        nElectronsFilter = (ak.count_nonzero(abs(events.lepPdgId)==electron_pdgId, axis=1)==1) & (events.nLeptons==1)
        cscClusterFilter = events.nCscRechitClusters==1
        if fromTau:
            fromTauFilter = ak.flatten(events.gTauEDecay)
        else:
            fromTauFilter = ak.ones_like(events.runNum, dtype=bool)
        event_mask = ak.zip({"noiseFilters":noiseFilters, "triggerFilter":triggerFilter, 
                             "nElectronsFilter":nElectronsFilter, "cscClusterFilter":cscClusterFilter, 
                             "fromTauFilter":fromTauFilter})
        if not mask_to_exclude==None:
            event_mask = ak.without_field(event_mask, mask_to_exclude)

        return event_mask
    
    def cscClusterSelections(self, cscClusters, mask_to_exclude=None):
        initial_mask = ak.ones_like(cscClusters.cscRechitClusterSize, dtype=bool)
        cluster_mask = ak.zip({"initial_mask":initial_mask, "muonVeto": cscClusters.cscRechitClusterMuonVetoPt<30, "jetVeto": cscClusters.cscRechitClusterJetVetoPt<10, "cluster_size": cscClusters.cscRechitClusterSize>=100, "inTime": (cscClusters.cscRechitClusterTimeWeighted>-5) & (cscClusters.cscRechitClusterTimeWeighted<12.5) & (cscClusters.cscRechitClusterNRechitME1112==0)})
        #cluster_mask = ak.zip({"initial_mask":initial_mask})
        if self.isMC and self.applyGenInfo:
            cluster_mask = ak.with_field(cluster_mask, cscClusters.cscRechitCluster_match_gLLP, "match_gLLP")
            cluster_mask = ak.with_field(cluster_mask, (abs(cscClusters.cscRechitCluster_match_gLLP_decay_z)>400) & (abs(cscClusters.cscRechitCluster_match_gLLP_decay_z)<1100), "LLP_inCSCs")
        
        if not mask_to_exclude==None:
            cluster_mask = ak.without_field(cluster_mask, mask_to_exclude)
        
        return cluster_mask

    def electronSelections(self, electrons, mask_to_exclude=None, invertElectronId=False):
        
        electron_ID_mask = electrons.lepTightId==(not invertElectronId) #change back when not looking at the cluster size!
        electron_ISO_mask = electrons.lepPassTightIso
        electron_pT_mask = electrons.lepPt>5
        electron_mask = ak.zip({"electron_ID_mask":electron_ID_mask, "electron_pT_mask":electron_pT_mask, "electron_ISO_mask":electron_ISO_mask})
        
        if not mask_to_exclude==None:
            electron_mask = ak.without_field(electron_mask, mask_to_exclude)

        return electron_mask

    def gLepSelections(self, gLeptons):
        return ak.zip({"no_cut": ak.ones_like(gLeptons.gLepPt, dtype=bool)})
    
    def gLLPSelections(self, gLLPs):
        return ak.zip({"no_cut": ak.ones_like(gLLPs.gLLP_pt, dtype=bool)})

    
    
    
    
    
    ###     PROCESSOR   ####
    ########################
    def process(self, events, hists_to_process: list=None, fillGenHists = True, fromTau = False, electronID: str='', branchNames_for_invertElectronID = ["CSC_Cluster_Size"]):
        
        hist_list = []
        if hists_to_process==None and self.isMC: #this means process all hists, since no specific ones are passed to the processor
            hist_list = [*list(self.eventLevel_hists_dict.keys()), *list(self.cscCluster_hists_dict.keys()), *list(self.gLLP_hists_dict.keys()), *list(self.gLepton_hists_dict.keys()), *list(self.electron_hists_dict.keys())]
        elif hists_to_process==None:
            hist_list = [*list(self.eventLevel_hists_dict.keys()), *list(self.cscCluster_hists_dict.keys()), *list(self.electron_hists_dict.keys())]
        else:
            hist_list = hists_to_process
        if not self.make_ElectronID_hists:
            hist_list = [hist for hist in hist_list if "ID" not in hist]
        #print(hist_list)
        output = self.accumulator.copy()
        
        #generate analysis objects
        electrons = self.buildRecoElectrons(events)
        cscClusters = self.buildCscRechitClusters(events)
        
        #apply cuts to analysis objects and do event-level selections
        cscCluster_mask = self.cscClusterSelections(cscClusters)
        electron_mask = self.electronSelections(electrons)
        electron_mask_invertId = self.electronSelections(electrons, invertElectronId=True)
        event_mask = self.eventSelections(events, fromTau=fromTau)

        #print(electrons['lepPt'].compute())

        if self.isMC and self.applyGenInfo:
            gLeptons = self.buildGLeptons(events)
            gLLP = self.buildGLLP(events)
            gLLP_mask = self.gLLPSelections(gLLP)
            gLepton_mask = self.gLepSelections(gLeptons)
            
            total_mask_electron = buildElectronMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask, gLLPMask=gLLP_mask, gLeptonMask=gLepton_mask, MC=True)
            #change back cluster mask
            total_mask_cscCluster = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask, gLLPMask=gLLP_mask,  gLeptonMask=gLepton_mask, MC=True)
            total_mask_cscCluster_invertElectronId = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask_invertId, gLLPMask=gLLP_mask,  gLeptonMask=gLepton_mask, MC=True)
            total_mask_event = buildEventMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask, gLLPMask=gLLP_mask,  gLeptonMask=gLepton_mask, MC=True)
            total_mask_gLLP = combineMasks(gLLP_mask)
            total_mask_gLepton = combineMasks(gLepton_mask)
            
        
        else:
            #use dummy true mask for gen-level masks
            total_mask_gLLP = None
            total_mask_gLepton = None
            total_mask_electron = buildElectronMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask, MC=False)
            total_mask_cscCluster_invertElectronId = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask_invertId, MC=False)
            total_mask_cscCluster = buildClusterMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask, MC=False)
            total_mask_event = buildEventMask(event_mask=event_mask, cluster_mask=cscCluster_mask, electron_mask=electron_mask, MC=False)

        #print(total_mask_event.compute())
        #print("Generated Masks, Starting to Fill Event-Level Histograms")

        collection_mask_mapping = {"electron": self.electronSelections, "cscRechitCluster": self.cscClusterSelections, "event": self.eventSelections}
        collection_objects = {"electron": electrons, "cscRechitCluster": cscClusters, "event": events}
        collection_totalMaskFunction_mapping = {"electron": buildElectronMask, "cscCluster": buildClusterMask, "event": buildEventMask}
        collection_mask_dict = {"electron": total_mask_electron, "cscCluster": total_mask_cscCluster, "event": total_mask_event,
                                "gLLP": total_mask_gLLP, "gLepton": total_mask_gLepton}

        def computeModMask(excluded_mask_collection, excluded_mask):
            mod_mask = collection_mask_mapping[excluded_mask_collection](collection_objects[excluded_mask_collection], excluded_mask)
            mask_dict_tmp = copy.deepcopy(collection_mask_dict[excluded_mask_collection])
            mask_dict_tmp[excluded_mask_collection] = mask_dict_tmp
            event_mask = collection_totalMaskFunction_mapping[excluded_mask_collection](**mask_dict_tmp)
            return event_mask


        to_fill = {}

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
            #eventLevelHist.fill(plot=events[eventLevelHist.metadata["branch"]][hist_mask])
            to_fill[plot] = (eventLevelHist, events[eventLevelHist.metadata["branch"]][hist_mask])
            #output[plot] = eventLevelHist


        if self.isMC and fillGenHists and self.applyGenInfo:
            #("filling gen-level hists")
            for plot in self.gLLP_hists_dict.keys():
                if plot not in hist_list:continue
                gLLPHist = self.gLLP_hists_dict[plot]
                added_mask = (events[gLLPHist.metadata["mask_branch"]]>=gLLPHist.metadata["mask_lowVal"]) & (events[gLLPHist.metadata["mask_branch"]]<=gLLPHist.metadata["mask_highVal"])
                if gLLPHist.metadata["mask_collection"]!="event":
                    added_mask = (events[collection_branch_mapping[gLLPHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
                hist_mask = (added_mask) & (total_mask_gLLP)
                #gLLPHist.fill(plot=ak.flatten(events[gLLPHist.metadata["branch"]][hist_mask]))
                to_fill[plot] = (gLLPHist, ak.flatten(events[gLLPHist.metadata["branch"]][hist_mask]))
                #output[plot] = gLLPHist
            
            for plot in self.gLepton_hists_dict.keys():
                if plot not in hist_list:continue
                gLepHist = self.gLepton_hists_dict[plot]
                added_mask = (events[gLepHist.metadata["mask_branch"]]>=gLepHist.metadata["mask_lowVal"]) & (events[gLepHist.metadata["mask_branch"]]<=gLepHist.metadata["mask_highVal"])
                if gLepHist.metadata["mask_collection"]!="event":
                    (events[collection_branch_mapping[gLepHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
                hist_mask = (added_mask) & (total_mask_gLepton)
                #gTauHist.fill(plot=ak.flatten(events[gTauHist.metadata["branch"]][hist_mask]))
                to_fill[plot] = (gLepHist, ak.flatten(events[gLepHist.metadata["branch"]][hist_mask]))
                #output[plot] = gTauHist
            
        
        #print("filling reco electron hists")
        for plot in self.electron_hists_dict.keys():
            if plot not in hist_list:continue
            electronHist = self.electron_hists_dict[plot]
            added_mask = (events[electronHist.metadata["mask_branch"]]>=electronHist.metadata["mask_lowVal"]) & (events[electronHist.metadata["mask_branch"]]<=electronHist.metadata["mask_highVal"])
            if electronHist.metadata["mask_collection"]!="event":
                added_mask = (events[collection_branch_mapping[electronHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
            hist_mask = (added_mask) & (total_mask_electron)
            #tauHist.fill(plot=ak.flatten(events[tauHist.metadata["branch"]][hist_mask]))
            to_fill[plot] = (electronHist, ak.flatten(events[electronHist.metadata["branch"]][hist_mask]))
            #output[plot] = tauHist
        
        for plot in self.cscCluster_hists_dict.keys():
            if plot not in hist_list:continue
            #print(plot)
            cscClusterHist = self.cscCluster_hists_dict[plot]
            added_mask = (events[cscClusterHist.metadata["mask_branch"]]>=cscClusterHist.metadata["mask_lowVal"]) & (events[cscClusterHist.metadata["mask_branch"]]<=cscClusterHist.metadata["mask_highVal"])
            if cscClusterHist.metadata["mask_collection"]!="event":
                added_mask = (events[collection_branch_mapping[cscClusterHist.metadata["mask_collection"]]]>0) & (combineMasksFlatten(ak.zip({"mask":added_mask})))
            hist_mask = (added_mask) & (total_mask_cscCluster)
            if cscClusterHist.metadata['Mask_To_Exclude']!=None:
                #print("excluding mask")
                cluster_mask = computeModMask(cscClusterHist.metadata['Mask_To_Exclude_Collection'], cscClusterHist.metadata['Mask_To_Exclude'])
            elif plot=="CSC_Cluster_Size" and not self.isMC:
                #print("Inverting Electron TightID in data for cluster size")
                cluster_mask = total_mask_cscCluster_invertElectronId
            else:
                cluster_mask = total_mask_cscCluster
            hist_mask = (added_mask) & (cluster_mask)
            to_fill[plot] = (cscClusterHist, ak.flatten(events[cscClusterHist.metadata["branch"]][hist_mask]))
            #output[plot] = cscClusterHist
        
        client = Client(memory_limit="12GB", n_workers=1, 
                threads_per_worker=1, local_directory="/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/dask_temp")
        #print("about to compute")
        
        computed_arrays = dask.compute(*[v[1] for v in to_fill.values()])

        for i, plot in enumerate(to_fill):
            hist_obj, _ = to_fill[plot]
            hist_obj.fill(plot=computed_arrays[i])
            output[plot] = hist_obj

        
        client.close()
        
        
        return output



    def postprocess(self, accumulator):
        return accumulator
