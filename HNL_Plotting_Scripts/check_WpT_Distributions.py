import sys
import matplotlib.pyplot as plt
sys.path.insert(0,"../python/HNL_Plotting_HelperFunctions")
import MuonSystemReader
import awkward as ak
import numpy as np
import scipy




new_MC= "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_021726/HNL_tau_mN_2_ctau_100/normalized/HNL_tau_mN_2_ctau_100_109080pb_weighted.root"
new_events  = MuonSystemReader.loadTree_nanoFactory(new_MC)
newPt = np.array(new_events.W_B_genPt.compute())

old_MC = "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_021726/HNL_tau_mN_2_ctau_1000/normalized/HNL_tau_mN_2_ctau_1000_109080pb_weighted.root"
old_events = MuonSystemReader.loadTree_nanoFactory(old_MC)
oldPt = np.array(old_events.W_B_genPt.compute())

import uproot
reference_MC = "root://cmseos.fnal.gov//store/group/lpclonglived/displacedJetMuonNtuple/MDSNano/Run3Summer24/v2/WtoENu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_WtoENu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8_v2/251204_161951/0000/EXO-RunIII2024Summer24NanoAODv15-00307_1.root"
with uproot.open(reference_MC) as file:
    # Syntax: file["TreeName"]["BranchName"]
    pt = file["Events"]["GenPart_pt"].array()
    
    pdgId = file["Events"]["GenPart_pdgId"].array()
    #print(pdgId[0])
    # Convert the branch data into a NumPy array

    pt = pt[np.abs(pdgId)==24]
    pt = pt[[len(sub) > 0 for sub in pt]]
    pt = pt[:,-1]
    print(pt)



bins = np.arange(0,200,5)
plt.hist(newPt, bins=bins, histtype='step', label="updated signal MC", density=True)
plt.hist(oldPt, bins=bins, histtype='step', label="old signal MC", density=True)
plt.hist(pt, bins=bins, histtype='step', label="WtoENu MC file", density=True)

plt.legend()
plt.savefig("pT_weighting/test.png")







