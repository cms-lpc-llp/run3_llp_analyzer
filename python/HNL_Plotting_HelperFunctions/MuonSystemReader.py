#!/usr/bin/env python3
import coffea
import awkward as ak
import dask
import dask_awkward as dak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
#import uproot
import uproot
uproot.source.xrootd.use_vector_read = False


import os


repoPath = os.environ['CMSSW_BASE']+'/src/run3_llp_analyzer/'

trigger_paths = {
    'electron': 'HLT_CscCluster100_Ele5',
    'muon': 'HLT_CscCluster100_Mu5',
    'tau': 'HLT_CscCluster100_PNetTauhPFJet10_Loose'
}


def reshape_DNN_Scores(DNN_Scores, nClusters):
    DNN_Scores_flat = ak.flatten(DNN_Scores)
    DNN_Scores_flat = DNN_Scores_flat[DNN_Scores_flat>=0]
    resize_DNN = ak.unflatten(DNN_Scores_flat, nClusters)
    return resize_DNN

def loadTree_nanoFactory(f, isMC=True, trigger="tau"):
    events = NanoEventsFactory.from_root(f,schemaclass=BaseSchema, metadata={"dataset":"HNL" }, steps_per_file=4).events()
    #if not isMC: 
    #    events = events[events[trigger_paths[trigger]]]
        #if trigger=='tau':
        #    events = events[(events['nTaus']>0) & (events['nCscRechitClusters']>0) & (ak.any(events['tauIsVLoose'], axis=1))]
    
    DNN_branch = events["cscRechitClusterDNN_bkgMC_plusBeamHalo"]
    events = ak.without_field(events, "cscRechitClusterDNN_bkgMC_plusBeamHalo")
    nClusters = events["nCscRechitClusters"]
    DNN_branch = reshape_DNN_Scores(DNN_branch, nClusters)
    
    events = ak.with_field(events, DNN_branch, "cscRechitClusterDNN_bkgMC_plusBeamHalo")#.persist()
    return events
    

def loadTree_dask(f, treePath="MuonSystem"):
    return uproot.dask(repoPath+f, library="ak")

