#!/bin/sh

tau_label="mN_{#tau}=2GeV ctau=1m"
e_label="mN_{e}=2GeV ctau=1m"
data_label="2024 Data"

CUDA_STUB_PATH="/cvmfs/cms.cern.ch/el8_amd64_gcc12/external/cuda/12.4.1-fc5cb0e72dba64b6abbf00089f3a044c/lib64/stubs"
if [ -n "${LD_LIBRARY_PATH}" ]; then
  new_ld=""
  old_ifs=$IFS
  IFS=:
  for p in $LD_LIBRARY_PATH; do
    if [ "$p" != "$CUDA_STUB_PATH" ]; then
      if [ -z "$new_ld" ]; then
        new_ld="$p"
      else
        new_ld="${new_ld}:$p"
      fi
    fi
  done
  IFS=$old_ifs
  export LD_LIBRARY_PATH="$new_ld"
fi

export PYTHONWARNINGS="ignore:.*coffea\\.nanoevents\\.methods\\.vector will be removed.*:FutureWarning,ignore:Please ensure that dask\\.awkward<.* is partitionwise-compatible.*:UserWarning"

# echo "cut based first"
# python3 run_HNL_Processor.py \
#   --flavor e \
#   --cut-based-id \
#   --apply-gen-info \
#   --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_012026/HNL_e_mN_2_ctau_1000/normalized/HNL_e_mN_2_ctau_1000_109080pb_weighted.root \
#   --label "${e_label}" \
#   --is-mc true \
#   --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_012026/HNL_tau_mN_2_ctau_1000/normalized/HNL_tau_mN_2_ctau_1000_109080pb_weighted.root \
#   --label "${tau_label}" \
#   --is-mc true \
#   --sample  root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/2024_Data_e_analyzer_update_012026/ \
#   --label "${data_label}" \
#   --is-mc false \
#    --hists electron_dz \
#   --outdir e_category_cut-based_MET_dz


echo "MVA based second"
python3 run_HNL_Processor.py \
  --flavor e \
  --no-cut-based-id \
  --apply-gen-info \
  --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_021726/HNL_e_mN_2_ctau_1000/normalized/HNL_e_mN_2_ctau_1000_109080pb_weighted.root \
  --label "${e_label}" \
  --is-mc true \
  --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/analyzer_update_021726/HNL_tau_mN_2_ctau_1000/normalized/HNL_tau_mN_2_ctau_1000_109080pb_weighted.root \
  --label "${tau_label}" \
  --is-mc true \
  --sample  root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/2024_Data_e_analyzer_update_021726/ \
  --label "${data_label}" \
  --is-mc false \
  --hists HT_SoftJets_LargerThan3p139 \
  --outdir e_category_MVA_PT_MET_ClusterSize160_IPCuts
