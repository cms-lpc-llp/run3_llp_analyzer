import sys
import matplotlib.pyplot as plt
sys.path.insert(0,"../python/HNL_Plotting_HelperFunctions")
import MuonSystemReader
import awkward as ak
import numpy as np
import scipy

def expo(x, c=100000, width=0.01):
    return c*np.exp(-width*x)

sample_ctau=100
reweight_ctau=1000


HNL_2GeV_10ctau_path = "root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_021726/HNL_tau_mN_2_ctau_100/normalized/HNL_tau_mN_2_ctau_100_109080pb_weighted.root"
signal_events  = MuonSystemReader.loadTree_nanoFactory(HNL_2GeV_10ctau_path)

ctaus = np.array(ak.flatten(signal_events.gLLP_ctau.compute()))

weights = sample_ctau/reweight_ctau*np.exp(10*ctaus*(1/sample_ctau-1/reweight_ctau))
#print(ak.sum(weights).compute())
weights = weights/ak.sum(weights)*len(ctaus)
print(max(weights))
print(np.sum(weights))
print(len(ctaus))

#ctaus = ctaus*weights

ctaus = np.array(ak.flatten(signal_events.gLLP_ctau.compute()))
bin_edges = np.arange(0, reweight_ctau, reweight_ctau/100)
counts,edges = np.histogram(ctaus, bin_edges, weights=weights)
centers = edges[:-1] + np.diff(edges) / 2

popt, pcov = scipy.optimize.curve_fit(expo, centers,counts, p0=[80000, 0.01])

plt.hist(ctaus, bins = bin_edges, weights=weights, histtype='step')
x = np.linspace(centers[0], centers[-1], 1000)
y = expo(x, *popt)
plt.xlim(0,reweight_ctau)
#plt.yscale('log')
#plt.xscale('log')
plt.plot(x,y, label=f"Fit: ctau = {1/popt[1]:.2f} cm")
plt.legend()
plt.savefig(f"check_lifetime_reweight_plots/ctaus_weighted{reweight_ctau}.png")
plt.clf()
plt.hist(weights, histtype='step')
plt.savefig(f"check_lifetime_reweight_plots/weights{reweight_ctau}.png")


