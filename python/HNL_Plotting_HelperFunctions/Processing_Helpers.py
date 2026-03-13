#!/usr/bin/env python3
import coffea
import awkward as ak
from coffea import processor

from coffea.nanoevents.methods import vector
from coffea.nanoevents import NanoEventsFactory, BaseSchema
import uproot
import numpy as np

import sys

import HNL_Processor_v2_mu
sys.path.append(".")
import HNL_Processor
import HNL_Processor_v2
import HNL_Processor_v2_e
import MuonSystemReader
import os

data_samples = [
    "Muon0-Run2024B-PromptReco-v1",
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

files_XROOTD_Error_Data = { 'e': [
    
],

'mu': [
   
],

'tau': []
}

path_to_bad_files = {'e': '/src/run3_llp_analyzer/HNL_Plotting_Scripts/local_copies/e/',
                     'mu': '/src/run3_llp_analyzer/HNL_Plotting_Scripts/local_copies/mu/',
                    'tau': '/src/run3_llp_analyzer/HNL_Plotting_Scripts/local_copies/tau/'}

def deltaPhi(tau_cscCluster_phi_diffs):
    '''
    deltaPhi returned where input array does not lie in (-pi, pi) range
    '''
    while(ak.count_nonzero(tau_cscCluster_phi_diffs>np.pi).compute()>0):
            tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs<np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs-(2*np.pi))
    while(ak.count_nonzero(tau_cscCluster_phi_diffs<-1*np.pi).compute()>0):
        tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs>-1*np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs+(2*np.pi))
    return tau_cscCluster_phi_diffs

def deltaPhiAk(tau_cscCluster_phi_diffs):
    '''
    deltaPhi returned where input array does not lie in (-pi, pi) range
    '''
    while(ak.count_nonzero(tau_cscCluster_phi_diffs>np.pi)>0):
            tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs<np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs-(2*np.pi))
    while(ak.count_nonzero(tau_cscCluster_phi_diffs<-1*np.pi)>0):
        tau_cscCluster_phi_diffs = ak.where(tau_cscCluster_phi_diffs>-1*np.pi, tau_cscCluster_phi_diffs, tau_cscCluster_phi_diffs+(2*np.pi))
    return tau_cscCluster_phi_diffs


def processData(data_path_base: str, hists_to_process: list=None, tau_cluster_topo_hists=True, trigger='tau'):

    data_events_list = [data_path_base+sample+"/normalized/"+sample+"_goodLumi.root" for sample in data_samples]
    
    outputs = []
    for sample in data_events_list:
        print(sample)
        data_events  = MuonSystemReader.loadTree_nanoFactory(sample, isMC=False, trigger=trigger)
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

def processData_v2(data_path_base: str, hists_to_process: list=None, trigger='tau', bin_multiplier=1):

    data_events_list = [data_path_base+sample+"/normalized/"+sample+"_goodLumi.root" for sample in data_samples]
    
    outputs = []
    for sample in data_events_list:
        correction_needed = False
        for file in files_XROOTD_Error_Data['tau']:
            if file in sample:
                correction_needed = True
                break
        if correction_needed:
            sample_corrected = sample.replace(data_path_base+sample+"/normalized/", os.environ['CMSSW_BASE']+path_to_bad_files['tau'])
        else:
            sample_corrected = sample
        print(sample_corrected)
        data_events  = MuonSystemReader.loadTree_nanoFactory(sample_corrected, isMC=False, trigger=trigger)
        processor_data = HNL_Processor_v2.HNL_Processor_v2(bin_multiplier = bin_multiplier)
        outputs.append(processor_data.process(data_events, hists_to_process = hists_to_process))

    output_data = {}
    for idx, output in enumerate(outputs):
        for key, hist in output.items():
            if idx==0:
                output_data[key] = hist
            else:
                output_data[key] = output_data[key]+hist

    return output_data

def processData_v2_e(data_path_base: str, hists_to_process: list=None, trigger='electron', bin_multiplier=1):

    data_events_list = [data_path_base+sample+"/normalized/"+sample+"_goodLumi.root" for sample in data_samples]
    
    outputs = []
    for sample in data_events_list:
        correction_needed = False
        for file in files_XROOTD_Error_Data['e']:
            if file in sample:
                correction_needed = True
                break
        if correction_needed:
            sample_corrected = sample.replace(data_path_base+sample+"/normalized/", os.environ['CMSSW_BASE']+path_to_bad_files['e'])
        else:
            sample_corrected = sample
        print(sample_corrected)
        data_events  = MuonSystemReader.loadTree_nanoFactory(sample_corrected, isMC=False, trigger=trigger)
        processor_data = HNL_Processor_v2_e.HNL_Processor_v2_e(bin_multiplier = bin_multiplier)
        outputs.append(processor_data.process(data_events, hists_to_process = hists_to_process))

    output_data = {}
    for idx, output in enumerate(outputs):
        for key, hist in output.items():
            if idx==0:
                output_data[key] = hist
            else:
                output_data[key] = output_data[key]+hist

    return output_data

def processData_v2_mu(data_path_base: str, hists_to_process: list=None, trigger='muon', bin_multiplier = 1):

    data_events_list = [data_path_base+sample+"/normalized/"+sample+"_goodLumi.root" for sample in data_samples]
    
    outputs = []
    for sample in data_events_list:
        correction_needed = False
        for file in files_XROOTD_Error_Data['mu']:
            if file in sample:
                correction_needed = True
                break
        if correction_needed:
            sample_corrected = sample.replace(data_path_base+sample+"/normalized/", os.environ['CMSSW_BASE']+path_to_bad_files['mu'])
        else:
            sample_corrected = sample
        print(sample_corrected)
        data_events  = MuonSystemReader.loadTree_nanoFactory(sample_corrected, isMC=False, trigger=trigger)
        processor_data = HNL_Processor_v2_mu.HNL_Processor_v2_mu(bin_multiplier = bin_multiplier)
        outputs.append(processor_data.process(data_events, hists_to_process = hists_to_process))

    output_data = {}
    for idx, output in enumerate(outputs):
        for key, hist in output.items():
            if idx==0:
                output_data[key] = hist
            else:
                output_data[key] = output_data[key]+hist

    return output_data

def make_data_list(data_path_base, flavor="tau", XROOTD_Replace=True):
    '''
    return lists with data files for each era to read from
    if XROOTD_Replace is set to true, then the path to local versions of problematic files to read through XROOTD will be pointed too
    '''
    data_events_list = [data_path_base+sample+"/normalized/"+sample+"_goodLumi.root" for sample in data_samples]
    if XROOTD_Replace:
        corrected_list = []
        for sample in data_events_list:
            correction_needed = False
            for file in files_XROOTD_Error_Data[flavor]:
                if file in sample:
                    correction_needed = True
                    break
            if correction_needed:
                sample_corrected = sample.replace(sample[:sample.rfind('/')], os.environ['CMSSW_BASE']+path_to_bad_files[flavor])
            else:
                sample_corrected = sample
            corrected_list.append(sample_corrected)
    else:
        return data_events_list

    return corrected_list


