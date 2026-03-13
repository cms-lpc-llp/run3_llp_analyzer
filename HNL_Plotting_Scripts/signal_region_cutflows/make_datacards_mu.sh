python3 Run_Signal_Ctau_Scan_ABCD.py \
  --cutflow mu_SR/cuts_sensitivity_mu_DNN_dz_dxy_noMET_MC.yaml \
  --sample signal_tau_2GeV_100mm --signal-name HNL_tau_mu_optimized --sample-ctau 100 \
  --sample-config sample_configs/mu_signal/HNL_mu_samples.yaml \
  --data-abcd SR_Cutflows/mu_SR__cuts_sensitivity_mu_DNN_dz_dxy_noMET/Data/abcd.txt \
  --ctaus 10,25,50,100,250,500,750,1000,2500,5000,7500,10000 \
  --mass 2 \
  --outdir SR_Cutflows \
  --datacards-outdir /uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit/Datacards/first_ABCD_counts/ \
  --size-cut 184
#--sample signal_mu_2GeV_1000mm --signal-name HNL_mu_optimzied --sample-ctau 1000 \
