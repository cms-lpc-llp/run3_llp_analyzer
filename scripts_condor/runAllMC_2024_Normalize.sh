#!/bin/bash

echo "Submitting jobs e-HNL signal"
python3 scripts_condor/submit_normalize_LPC.py HNL_e_mN_2_ctau_1000 $1 e-HNL_$1

echo "Submitting jobs mu-HNL signal"
python3 scripts_condor/submit_normalize_LPC.py HNL_mu_mN_2_ctau_1000 $1 mu-HNL_$1

echo "Submitting jobs tau-HNL signal"
python3 scripts_condor/submit_normalize_LPC.py HNL_tau_mN_2_ctau_1000 $1 tau-HNL_$1

echo "Submitting jobs WtoENu MC Bkg"
python3 scripts_condor/submit_normalize_LPC.py WtoENu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 WtoENu_$1

echo "Submitting jobs WtoMuNu MC Bkg"
python3 scripts_condor/submit_normalize_LPC.py WtoMuNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 WtoMuNu_$1

echo "Submitting jobs WtoTauNu MC Bkg"
python3 scripts_condor/submit_normalize_LPC.py WtoTauNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 WtoTauNu_$1

###TODO: Z->ee bkg


echo "Submitting jobs ZtoMuMu MC Bkg"
python3 scripts_condor/submit_normalize_LPC.py DYto2Mu_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8 $1 ZtoMuMu_$1

echo "Submitting jobs ZtoTauTau MC Bkg"
python3 scripts_condor/submit_normalize_LPC.py DYto2Tau-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 ZtoTauTau_$1
