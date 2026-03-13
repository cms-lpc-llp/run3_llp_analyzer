#!/bin/bash

# echo "Submitting jobs e-HNL signal, ctau 1000 mN 2"
# python3 scripts_condor/submit_condor_LPC.py HNL_e_mN_2_ctau_1000 mdsnano_hnl 1 $1 e-HNL_$1 10 4096

# echo "Submitting jobs mu-HNL signal, ctau 1000 mN 2"
# python3 scripts_condor/submit_condor_LPC.py HNL_mu_mN_2_ctau_1000 mdsnano_hnl 2 $1 mu-HNL_$1 10 4096

# echo "Submitting jobs tau-HNL signal, ctau 1000 mN 2"
# python3 scripts_condor/submit_condor_LPC.py HNL_tau_mN_2_ctau_1000 mdsnano_hnl 3 $1 tau-HNL_$1 10 4096

# echo "Submitting jobs tau-HNL signal, ctau 100 mN 2"
# python3 scripts_condor/submit_condor_LPC.py HNL_tau_mN_2_ctau_100 mdsnano_hnl 3 $1 tau-HNL_100_ctau_$1 2 4096

# echo "Submitting jobs tau-HNL signal, ctau 100 mN 1"
# python3 scripts_condor/submit_condor_LPC.py HNL_tau_mN_1_ctau_100 mdsnano_hnl 3 $1 tau-HNL_$1 2 4096

# echo "Submitting jobs tau-HNL signal, ctau 100 mN 3"
# python3 scripts_condor/submit_condor_LPC.py HNL_tau_mN_3_ctau_100 mdsnano_hnl 3 $1 tau-HNL_$1 2 4096


echo "Submitting jobs tau-HNL signal from Bs, ctau 1000 mN 2"
python3 scripts_condor/submit_condor_LPC.py HNL_B_tau_mN_2_ctau_1000 mdsnano_hnl 13 $1 tau-HNL_fromB_2GeV_ctau1000_try2_$1 2 4096


# echo "Submitting jobs tau-HNL signal from Bs, ctau 1000 mN 3"
# python3 scripts_condor/submit_condor_LPC.py HNL_B_tau_mN_3_ctau_1000 mdsnano_hnl 13 $1 tau-HNL_fromB_3GeV_ctau1000_try2_$1 2 4096

# echo "Submitting jobs tau-HNL signal from Bs, ctau 10000 mN 3"
# python3 scripts_condor/submit_condor_LPC.py HNL_B_tau_mN_3_ctau_10000 mdsnano_hnl 13 $1 tau-HNL_fromB_3GeV_ctau10000_try2_$1 2 4096

# echo "Submitting jobs WtoENu MC Bkg"
# python3 scripts_condor/submit_condor_LPC.py WtoENu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 mdsnano_hnl 4 $1 WtoENu_$1 5 8000

# echo "Submitting jobs WtoMuNu MC Bkg"
# python3 scripts_condor/submit_condor_LPC.py WtoMuNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 mdsnano_hnl 5 $1 WtoMuNu_$1 5 8000

# echo "Submitting jobs WtoTauNu MC Bkg"
# python3 scripts_condor/submit_condor_LPC.py WtoTauNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8 mdsnano_hnl 6 $1 WtoTauNu_$1 5 8000

# ###TODO: Z->ee bkg


# echo "Submitting jobs ZtoMuMu MC Bkg"
# python3 scripts_condor/submit_condor_LPC.py DYto2Mu_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8 mdsnano_hnl 8 $1 ZtoMuMu_$1 5 8000

# echo "Submitting jobs ZtoTauTau MC Bkg"
# python3 scripts_condor/submit_condor_LPC.py DYto2Tau-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8 mdsnano_hnl 9 $1 ZtoTauTau_$1 5 8000
