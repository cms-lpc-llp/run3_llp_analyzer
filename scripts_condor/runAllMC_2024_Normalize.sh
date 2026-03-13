#!/bin/bash

#first argument: subdirectory after HNL_Tau_Search where the data sits
#second argument: integrated lumi to normalize to

# echo "Submitting jobs e-HNL signal"
# python3 scripts_condor/submit_normalize_LPC.py HNL_e_mN_2_ctau_1000 $1 e-HNL_$1 $2

# echo "Submitting jobs mu-HNL signal"
# python3 scripts_condor/submit_normalize_LPC.py HNL_mu_mN_2_ctau_1000 $1 mu-HNL_$1 $2

# echo "Submitting jobs tau-HNL signal"
# python3 scripts_condor/submit_normalize_LPC.py HNL_tau_mN_2_ctau_1000 $1 tau-HNL_$1 $2

# echo "Submitting jobs WtoENu MC Bkg"
# python3 scripts_condor/submit_normalize_LPC.py WtoENu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 WtoENu_$1 $2

# echo "Submitting jobs WtoMuNu MC Bkg"
# python3 scripts_condor/submit_normalize_LPC.py WtoMuNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 WtoMuNu_$1 $2

# echo "Submitting jobs WtoTauNu MC Bkg"
# python3 scripts_condor/submit_normalize_LPC.py WtoTauNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 WtoTauNu_$1 $2

# echo "Submitting jobs tau-HNL signal, ctau 100 mN 2"
# python3 scripts_condor/submit_normalize_LPC.py HNL_tau_mN_2_ctau_100 $1 tau-HNL_100_ctau_$1 $2

# echo "Submitting jobs tau-HNL signal, ctau 100 mN 1"
# python3 scripts_condor/submit_normalize_LPC.py HNL_tau_mN_1_ctau_100 $1 tau-HNL_100_1_$1 $2

# echo "Submitting jobs tau-HNL signal, ctau 100 mN 3"
# python3 scripts_condor/submit_normalize_LPC.py HNL_tau_mN_3_ctau_100 mdsnano_hnl $1 tau-HNL_100_3_$1 $2

echo "Submitting jobs tau-HNL signal from Bs, ctau 1000 mN 2"
python3 scripts_condor/submit_normalize_LPC.py HNL_B_tau_mN_2_ctau_1000 $1 tau-HNL_fromB_2GeV_ctau1000_$1 $2


# echo "Submitting jobs tau-HNL signal from Bs, ctau 1000 mN 3"
# python3 scripts_condor/submit_normalize_LPC.py HNL_B_tau_mN_3_ctau_1000 $1 tau-HNL_fromB_3GeV_ctau1000_$1 $2

# echo "Submitting jobs tau-HNL signal from Bs, ctau 10000 mN 3"
# python3 scripts_condor/submit_normalize_LPC.py HNL_B_tau_mN_3_ctau_10000 $1 tau-HNL_fromB_3GeV_ctau10000_$1 $2

###TODO: Z->ee bkg


# echo "Submitting jobs ZtoMuMu MC Bkg"
# python3 scripts_condor/submit_normalize_LPC.py DYto2Mu_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8 $1 ZtoMuMu_$1 $2

# echo "Submitting jobs ZtoTauTau MC Bkg"
# python3 scripts_condor/submit_normalize_LPC.py DYto2Tau-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8 $1 ZtoTauTau_$1 $2
