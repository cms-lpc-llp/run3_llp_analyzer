#!/bin/bash

# samples_list=(
#     Muon0-Run2024B-PromptReco-v1
#     Muon0-Run2024C-PromptReco-v1
#     Muon0-Run2024D-PromptReco-v1
#     Muon0-Run2024E-PromptReco-v1
#     Muon0-Run2024E-PromptReco-v2
#     Muon0-Run2024F-PromptReco-v1
#     Muon0-Run2024G-PromptReco-v1
#     Muon0-Run2024H-PromptReco-v1
#     Muon0-Run2024I-PromptReco-v1
#     Muon0-Run2024I-PromptReco-v2
#     Muon1-Run2024B-PromptReco-v1
#     Muon1-Run2024C-PromptReco-v1
#     Muon1-Run2024D-PromptReco-v1
#     Muon1-Run2024E-PromptReco-v1
#     Muon1-Run2024E-PromptReco-v2
#     Muon1-Run2024F-PromptReco-v1
#     Muon1-Run2024G-PromptReco-v1
#     Muon1-Run2024H-PromptReco-v1
#     Muon1-Run2024I-PromptReco-v1
#     Muon1-Run2024I-PromptReco-v2
# )

samples_list=(
    Muon1-Run2024H-PromptReco-v1
)

for sample in "${samples_list[@]}"; do
    echo "Submitting job for sample: $sample"
    python3 scripts_condor/submit_normalize_LPC.py "$sample" 2024_Data_e_retry 2024_Data_e_retry_normalize
done
echo "All jobs submitted for samples in the list."