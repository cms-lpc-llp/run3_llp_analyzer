#!/usr/bin/env python3

from pathlib import Path

import ABCD_VR1_ClosureStudies_Summary as base
import Processing_Helpers
import Produce_Cutflow_v2_forBs


_original_parse_args = base.parse_args


def parse_args():
    args = _original_parse_args()
    if args.outdir == "ABCD_VR1_ClosureStudies_Output":
        args.outdir = "ABCD_VR1_ClosureStudies_Output_forBs"
    return args


def _is_data_base_path(path: str) -> bool:
    if path.endswith("/"):
        return True
    lowered = path.lower()
    if lowered.endswith(".txt"):
        return False
    if lowered.startswith(("root://", "gsiftp://", "xrootd://")) and not lowered.endswith(".root"):
        return True
    return not lowered.endswith(".root")


def resolve_files_forbs(sample_cfg):
    if "data_path_base" in sample_cfg:
        flavor_key = sample_cfg.get("data_flavor") or base.flavor_to_data_key(sample_cfg.get("flavor", "tau"))
        return Processing_Helpers.make_data_list(sample_cfg["data_path_base"], flavor_key)

    files = sample_cfg.get("files") or sample_cfg.get("file_list") or sample_cfg.get("list_file")
    if not files:
        raise ValueError("Sample config must include 'files' or 'data_path_base'.")

    is_mc = base._to_bool(sample_cfg.get("isMC", sample_cfg.get("is_mc", False)), default=False)
    flavor = sample_cfg.get("flavor", "tau")

    if isinstance(files, str):
        if files.startswith("@"):
            return base.read_file_list(files[1:])
        if files.endswith(".txt") and Path(files).exists():
            return base.read_file_list(files)
        if not is_mc and _is_data_base_path(files):
            flavor_key = sample_cfg.get("data_flavor") or base.flavor_to_data_key(flavor)
            return Processing_Helpers.make_data_list(files, flavor_key)
        return [files]

    resolved = []
    for entry in files:
        if isinstance(entry, str) and entry.endswith(".txt") and Path(entry).exists():
            resolved.extend(base.read_file_list(entry))
        elif isinstance(entry, str) and not is_mc and _is_data_base_path(entry):
            flavor_key = sample_cfg.get("data_flavor") or base.flavor_to_data_key(flavor)
            resolved.extend(Processing_Helpers.make_data_list(entry, flavor_key))
        else:
            resolved.append(entry)
    return resolved


base.Produce_Cutflow_v2 = Produce_Cutflow_v2_forBs
base.parse_args = parse_args
base.resolve_files = resolve_files_forbs


if __name__ == "__main__":
    base.main()
