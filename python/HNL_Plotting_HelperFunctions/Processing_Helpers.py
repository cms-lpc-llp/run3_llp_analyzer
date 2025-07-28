#!/usr/bin/env python3
import coffea
import awkward as ak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot
import numpy as np

import sys
sys.path.append(".")
import HNL_Processor
import MuonSystemReader

data_samples = ["Muon0-Run2024B-PromptReco-v1",
    "Muon0-Run2024C-PromptReco-v1",
    "Muon0-Run2024D-PromptReco-v1",
    "Muon0-Run2024E-PromptReco-v1",
    "Muon0-Run2024E-PromptReco-v2",
    "Muon0-Run2024F-PromptReco-v1",
    "Muon0-Run2024G-PromptReco-v1",
    "Muon0-Run2024H-PromptReco-v1",
    "Muon0-Run2024I-PromptReco-v1",
    "Muon0-Run2024I-PromptReco-v2",
    "Muon1-Run2024B-PromptReco-v1",
    "Muon1-Run2024C-PromptReco-v1",
    "Muon1-Run2024D-PromptReco-v1",
    "Muon1-Run2024E-PromptReco-v1",
    "Muon1-Run2024E-PromptReco-v2",
    "Muon1-Run2024F-PromptReco-v1",
    "Muon1-Run2024G-PromptReco-v1",
    "Muon1-Run2024H-PromptReco-v1",
    "Muon1-Run2024I-PromptReco-v1",
    "Muon1-Run2024I-PromptReco-v2"
    ]


def deltaPhi(tau_cscCluster_phi_diffs):
    '''
    deltaPhi returned where input array does not like in (-pi, pi) range
    '''
    while(ak.count_nonzero(tau_cscCluster_phi_diffs>np.pi).compute()>0):
            tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs<np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs-(2*np.pi))
    while(ak.count_nonzero(tau_cscCluster_phi_diffs<-1*np.pi).compute()>0):
        tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs>-1*np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs+(2*np.pi))
    return tau_cscCluster_phi_diffs

def deltaPhiAk(tau_cscCluster_phi_diffs):
    '''
    deltaPhi returned where input array does not like in (-pi, pi) range
    '''
    while(ak.count_nonzero(tau_cscCluster_phi_diffs>np.pi)>0):
            tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs<np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs-(2*np.pi))
    while(ak.count_nonzero(tau_cscCluster_phi_diffs<-1*np.pi)>0):
        tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs>-1*np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs+(2*np.pi))
    return tau_cscCluster_phi_diffs


def processData(data_path_base: str, hists_to_process: list=None, tau_cluster_topo_hists=True):

    data_events_list = [data_path_base+sample+"/normalized/"+sample+"_goodLumi.root" for sample in data_samples]
    
    outputs = []
    for sample in data_events_list:
        print(sample)
        data_events  = MuonSystemReader.loadTree_nanoFactory(sample)
        processor_data = HNL_Processor.HNL_Processor()
        outputs.append(processor_data.process(data_events, hists_to_process = hists_to_process, tau_cluster_topo_hists=tau_cluster_topo_hists))

    output_data = {}
    for idx, output in enumerate(outputs):
        for key, hist in output.items():
            if idx==0:
                output_data[key] = hist
            else:
                output_data[key] = output_data[key]+hist

    return output_data
