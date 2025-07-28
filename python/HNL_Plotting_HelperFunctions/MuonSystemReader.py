#!/usr/bin/env python3
import coffea
import awkward as ak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot


repoPath = "/uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/run3_llp_analyzer/"




def loadTree_nanoFactory(f, isMC=True):
    return NanoEventsFactory.from_root(f,schemaclass=BaseSchema, metadata={"dataset":"HNL" }).events()
    
def loadTree_dask(f, treePath="MuonSystem"):
    return uproot.dask(repoPath+f, library="ak")

