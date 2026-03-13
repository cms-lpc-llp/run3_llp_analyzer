python3 Run_Signal_Ctau_Scan_ABCD.py \
  --cutflow tau_SR/cuts_sensitivity_tau_DNN_VVTight_MC.yaml \
  --sample signal_tau_2GeV_100mm --signal-name HNL_tauh_optimized --sample-ctau 100 \
  --sample-config sample_configs/tau_signal/HNL_tau_samples.yaml \
  --data-abcd SR_Cutflows/tau_SR__cuts_sensitivity_tau_DNN_VVTight/Data/abcd.txt \
  --ctaus 10,25,50,100,250,500,750,1000,2500,5000,7500,10000 \
  --mass 2 \
  --outdir SR_Cutflows \
  --datacards-outdir /uscms/home/amalbert/nobackup/el9_work/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit/Datacards/first_ABCD_counts/ \
  --size-cut 184
