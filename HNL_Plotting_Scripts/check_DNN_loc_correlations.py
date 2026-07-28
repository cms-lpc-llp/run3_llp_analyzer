import sys
import matplotlib.pyplot as plt
sys.path.insert(0,"../python/HNL_Plotting_HelperFunctions")
import MuonSystemReader
import awkward as ak
import numpy as np
import scipy




signal = "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_021726/HNL_tau_mN_2_ctau_100/normalized/HNL_tau_mN_2_ctau_100_109080pb_weighted.root"
signal_events  = MuonSystemReader.loadTree_nanoFactory(signal)
signal_events = signal_events[(signal_events.nCscRechitClusters==1)&(ak.any(signal_events.gLLP_csc, axis=1))]
signalDNN = np.array(ak.flatten(signal_events.cscRechitClusterDNN_bkgMC_plusBeamHalo).compute())
signalDNN = signalDNN[signalDNN>0]
clusterZ = np.array(ak.flatten(signal_events.cscRechitClusterZ).compute())
DNN_ranges = [(0.9,0.96),(0.96,0.99),(0.99,0.999),(0.999,0.99925),(0.99925,0.99975),(0.99975,1)]
bins = np.arange(500,1100,20 )
for DNN_range in DNN_ranges:
    clusterZ_masked = np.abs(clusterZ[np.logical_and(signalDNN>DNN_range[0],signalDNN<DNN_range[1])])
    plt.hist(clusterZ_masked,  bins=bins,histtype='step', label=str(DNN_range), density=True)
#plt.hist(oldPt, bins=bins, histtype='step', label="old signal MC", density=True)
#plt.hist(pt, bins=bins, histtype='step', label="WtoENu MC file", density=True)
plt.xlabel("|Cluster Z|")
plt.ylabel("Density")
plt.title("Cluster Z Distribution for DNN Slices in Signal")
plt.legend()
plt.savefig("DNN_study_clusterZ/clusterZDist_DNNSlices.png")







