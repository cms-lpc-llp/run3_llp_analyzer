## Environment Setup

Run `source setup.sh` to get `python 3.9`

## Get Cache of Ntuple and NanoAOD

Run `python3 wrapper_cache.py`, make sure to `base_path` and `list_path` point to the correct git repo path and list path for both ntuples and nanoAOD
Output cache are automatically saved to `run3_llp_analyzer/data/cache/`

## Produce Merging Commands

Run `python3 prepare_inputs.py` for data or `python3 prepare_inputs_mc.py` for MC which go through the cache for ntuple and nanoAOD to find files with overlapping events and generate text files in directory `merge_ntuples_commands` that can be used directly in condor jobs by `./MergeNtuples`

## Run MergeNtuples

Run `MergeNtuples` on condor which take the list of text files from the previous steps

```bash
cmsenv
python3 submit_MergeNtuple_caltech_cache.py
```
