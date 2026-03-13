#!/bin/sh

tau_label="mN_{#tau}=2GeV ctau=1m"
mu_label="mN_{#mu}=2GeV ctau=1m"
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

export PYTHONWARNINGS="ignore::FutureWarning:coffea\\.nanoevents\\.schemas\\..*,ignore::UserWarning:dask_awkward\\.lib\\.structure,ignore:Port 8787 is already in use.*:UserWarning:distributed\\.node,ignore:Creating scratch directories is taking a surprisingly long time.*:UserWarning"

echo "mu selection plots"
python3 run_HNL_Processor.py \
  --flavor mu \
  --not_normalized \
  --linscale \
  --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/2024_Data_mu_analyzer_update_012026/ \
  --label "${data_label}_looseIso" \
  --is-mc false \
  --iso loose \
  --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/2024_Data_mu_analyzer_update_012026/ \
  --label "${data_label}_tightIso" \
  --is-mc false \
  --iso tight \
  --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/2024_Data_mu_analyzer_update_012026/ \
  --label "${data_label}_vtightIso" \
  --is-mc false \
  --iso vtight \
  --sample root://cmseos.fnal.gov//store/group/lpclonglived/amalbert/HNL_Tau_Search/2024_Data_mu_analyzer_update_012026/ \
  --label "${data_label}_vvtightIso" \
  --is-mc false \
  --iso vvtight \
  --hists muon_pt \
  --outdir mu_category_MET_scanISO_zoomIn
  #--fine-only \
