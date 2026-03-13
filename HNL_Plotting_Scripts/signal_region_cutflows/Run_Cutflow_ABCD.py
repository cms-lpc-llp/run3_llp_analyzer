#!/usr/bin/env python3

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
import awkward as ak

repo_root = Path(__file__).resolve().parents[2]
helper_dir = repo_root / "python" / "HNL_Plotting_HelperFunctions"
sys.path.insert(0, str(repo_root))
sys.path.insert(0, str(helper_dir))

warnings.filterwarnings(
    "ignore",
    message=".*coffea\\.nanoevents\\.methods\\.vector will be removed.*",
    category=FutureWarning,
    module="coffea\\.nanoevents\\.schemas\\.fcc",
)
warnings.filterwarnings(
    "ignore",
    message="The value of the smallest subnormal for .*",
    category=UserWarning,
    module="numpy\\.core\\.getlimits",
)
warnings.filterwarnings(
    "ignore",
    message="You have set steps_per_file to 4,.*",
    category=RuntimeWarning,
    module="coffea\\.nanoevents\\.factory",
)
warnings.filterwarnings(
    "ignore",
    message="Please ensure that dask\\.awkward<.* is partitionwise-compatible.*",
    category=UserWarning,
    module="dask_awkward\\.lib\\.structure",
)

import MuonSystemReader
import Produce_Cutflow_v2
import Processing_Helpers
import analysis_helpers


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run signal region cutflows (i.e not closure studies) "
            "saving cutflows, ABCD tables"
        )
    )
    parser.add_argument(
        "--cutflow",
        required=True,
        help=(
            "Cutflow YAML config path relative to cuts_config (e.g. "
            "e_SR/cuts_sensitivity_study_DNN.yaml)."
        ),
    )
    parser.add_argument(
        "--sample",
        action="append",
        help=(
            "Sample name (if --sample-config is used) or an inline sample spec "
            "(key=value;key=value...)."
        ),
    )
    parser.add_argument(
        "--sample-config",
        default=None,
        help="YAML config that defines sample metadata keyed by name.",
    )
    parser.add_argument(
        "--size-cut",
        type=float,
        default=160.0,
        help="Baseline cluster size cut for ABCD (default: 160).",
    )
    parser.add_argument(
        "--dphi-cut",
        type=float,
        default=2.0,
        help="Baseline dPhi cut for ABCD (default: 2.0).",
    )
    parser.add_argument(
        "--blind",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Blind the ABCD signal bin (default: unblinded).",
    )
    parser.add_argument(
        "--outdir",
        default="SR_Cutflows",
        help="Output directory for results.",
    )
    parser.add_argument(
        "--sample-ctau",
        type=float,
        default=None,
        help="Override sample ctau used for MC lifetime reweighting.",
    )
    parser.add_argument(
        "--reweight-ctau",
        type=float,
        default=None,
        help="Override target ctau used for MC lifetime reweighting.",
    )
    return parser.parse_args()


def parse_list_arg(value: Optional[str], default: List[float]) -> List[float]:
    if not value:
        return list(default)
    return [float(x.strip()) for x in value.split(",") if x.strip()]


def read_file_list(path: str) -> List[str]:
    entries = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            entries.append(line)
    return entries


def normalize_flavor(flavor: str) -> str:
    norm = flavor.strip().lower()
    if norm in {"tau", "tauh", "had"}:
        return "Tau"
    if norm in {"mu", "muon", "m"}:
        return "Mu"
    if norm in {"e", "ele", "electron"}:
        return "Ele"
    return flavor


def flavor_to_trigger(flavor: str) -> str:
    norm = flavor.strip().lower()
    if norm in {"ele", "electron", "e"}:
        return "electron"
    if norm in {"mu", "muon", "m"}:
        return "muon"
    return "tau"


def flavor_to_data_key(flavor: str) -> str:
    norm = flavor.strip().lower()
    if norm in {"ele", "electron", "e"}:
        return "e"
    if norm in {"mu", "muon", "m"}:
        return "mu"
    return "tau"


def parse_inline_sample(spec: str) -> Dict[str, str]:
    if "=" not in spec:
        raise ValueError(
            "Inline sample specs must be key=value pairs separated by ';'. "
            "Example: name=signal_tau;files=lists/fewHNLFiles.txt;flavor=Tau;isMC=true"
        )
    out: Dict[str, str] = {}
    for part in spec.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"Invalid sample spec segment: {part}")
        key, value = part.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def load_sample_config(path: str) -> Dict[str, Dict]:
    with open(path, "r") as handle:
        data = yaml.safe_load(handle) or {}
    if "samples" in data:
        return data["samples"]
    return data


def _is_data_base_path(path: str) -> bool:
    if path.endswith("/"):
        return True
    lowered = path.lower()
    if lowered.endswith(".txt"):
        return False
    if lowered.startswith(("root://", "gsiftp://", "xrootd://")) and not lowered.endswith(".root"):
        return True
    return not lowered.endswith(".root")


def resolve_files(sample_cfg: Dict, *, is_mc: bool, flavor: str) -> List[str]:
    if "data_path_base" in sample_cfg:
        flavor_key = sample_cfg.get("data_flavor") or flavor_to_data_key(sample_cfg.get("flavor", "tau"))
        return Processing_Helpers.make_data_list(sample_cfg["data_path_base"], flavor_key)

    files = sample_cfg.get("files") or sample_cfg.get("file_list") or sample_cfg.get("list_file")
    if not files:
        raise ValueError("Sample config must include 'files' or 'data_path_base'.")

    if isinstance(files, str):
        if files.startswith("@"):
            return read_file_list(files[1:])
        if files.endswith(".txt") and os.path.exists(files):
            return read_file_list(files)
        if not is_mc and _is_data_base_path(files):
            flavor_key = sample_cfg.get("data_flavor") or flavor_to_data_key(flavor)
            return Processing_Helpers.make_data_list(files, flavor_key)
        return [files]

    resolved: List[str] = []
    for entry in files:
        if isinstance(entry, str) and entry.endswith(".txt") and os.path.exists(entry):
            resolved.extend(read_file_list(entry))
        elif isinstance(entry, str) and not is_mc and _is_data_base_path(entry):
            flavor_key = sample_cfg.get("data_flavor") or flavor_to_data_key(flavor)
            resolved.extend(Processing_Helpers.make_data_list(entry, flavor_key))
        else:
            resolved.append(entry)
    return resolved


def _flatten_branch(arr) -> np.ndarray:
    if arr is None:
        return np.array([])
    if hasattr(arr, "compute"):
        arr = arr.compute()
    try:
        return np.asarray(ak.flatten(arr))
    except Exception:
        return np.asarray(arr)


def _abcd_counts_from_branches(
    cluster_sizes: np.ndarray,
    dphi_vals: np.ndarray,
    size_cut: float,
    dphi_cut: float,
    normalization_factor: float,
    blind: bool,
) -> Tuple[float, float, Optional[float], Optional[float]]:
    mask_a = (cluster_sizes < size_cut) & (np.abs(dphi_vals) < dphi_cut)
    mask_b = (cluster_sizes < size_cut) & (np.abs(dphi_vals) >= dphi_cut)
    mask_c = (cluster_sizes >= size_cut) & (np.abs(dphi_vals) < dphi_cut)
    mask_d = (cluster_sizes >= size_cut) & (np.abs(dphi_vals) >= dphi_cut)

    A = float(np.count_nonzero(mask_a)) * normalization_factor
    B = float(np.count_nonzero(mask_b)) * normalization_factor
    C = float(np.count_nonzero(mask_c)) * normalization_factor
    D = float(np.count_nonzero(mask_d)) * normalization_factor

    A_unc = np.sqrt(A)
    B_unc = np.sqrt(B)
    C_unc = np.sqrt(C)
    D_unc = np.sqrt(D)

    if A <= 0 or B <= 0 or C <= 0:
        D_exp, D_exp_unc = np.nan, np.nan
    else:
        D_exp, D_exp_unc = analysis_helpers.compute_expected_signal_bin(A, A_unc, B, B_unc, C, C_unc)

    if blind:
        return D_exp, D_exp_unc, np.nan, np.nan
    return D_exp, D_exp_unc, D, D_unc



def cutflow_tag_from_path(cutflow_cfg: str) -> str:
    tag = cutflow_cfg.strip().replace("\\", "/")
    if tag.endswith(".yaml"):
        tag = tag[:-5]
    tag = tag.replace("/", "__")
    return tag or "cutflow"


def build_sample_defs(args: argparse.Namespace) -> Dict[str, Dict]:
    if args.sample_config:
        sample_map = load_sample_config(args.sample_config)
        resolved = {}
        names = args.sample or list(sample_map.keys())
        for name in names:
            if name not in sample_map:
                raise ValueError(f"Sample '{name}' not found in {args.sample_config}.")
            cfg = dict(sample_map[name])
            cfg["name"] = name
            resolved[name] = cfg
        return resolved

    if not args.sample:
        raise ValueError("--sample is required when --sample-config is not provided.")

    resolved = {}
    for spec in args.sample:
        cfg = parse_inline_sample(spec)
        name = cfg.get("name") or cfg.get("label")
        if not name:
            raise ValueError("Inline sample spec must include 'name'.")
        resolved[name] = cfg
    return resolved


def run_sample(
    name: str,
    cfg: Dict,
    cutflow_cfg: str,
    size_cut: float,
    dphi_cut: float,
    blind: bool,
    outdir: Path,
    sample_ctau_override: Optional[float] = None,
    reweight_ctau_override: Optional[float] = None,
) -> None:
    flavor = normalize_flavor(cfg.get("flavor", "Tau"))
    trigger = cfg.get("trigger") or flavor_to_trigger(flavor)
    is_mc = str(cfg.get("isMC", cfg.get("is_mc", False))).lower() in {"1", "true", "yes", "y", "t"}
    is_bkg = str(cfg.get("isBkg", cfg.get("is_bkg", False))).lower() in {"1", "true", "yes", "y", "t"}
    files = resolve_files(cfg, is_mc=is_mc, flavor=flavor)
    events = MuonSystemReader.loadTree_nanoFactory(files, isMC=is_mc, trigger=trigger)

    normalization_factor = 1.0
    if is_mc or float(cfg.get("signal_xsec", 1)) != 1:
        generated_events = ak.count(events.evtNum)
        gen_events_cfg = float(cfg.get("gen_events", 1))
        if gen_events_cfg == 1:
            num_events = generated_events.compute() if hasattr(generated_events, "compute") else float(generated_events)
        else:
            num_events = gen_events_cfg
        normalization_factor = Produce_Cutflow_v2.processed_lumi * float(cfg.get("signal_xsec", 1)) / num_events

    
    sample_ctau = (
        float(sample_ctau_override)
        if sample_ctau_override is not None
        else float(cfg.get("sample_ctau", 1000))
    )
    reweight_ctau = (
        float(reweight_ctau_override)
        if reweight_ctau_override is not None
        else float(cfg.get("reweight_ctau", 1000))
    )

    cutflow_df, abcd_df = Produce_Cutflow_v2.makeCutflow(
        events,
        cfg_file=cutflow_cfg,
        isMC=is_mc,
        ABCD=True,
        sizeCut=size_cut,
        dPhiCut=dphi_cut,
        blind=blind,
        flavor=flavor,
        sample_ctau=sample_ctau,
        reweight_ctau=reweight_ctau,
        signal_xsec=float(cfg.get("signal_xsec", 1)),
        genEvents=float(cfg.get("gen_events", 1)),
        isBkg = is_bkg,
    )

    cutflow_tag = cutflow_tag_from_path(cutflow_cfg)
    sample_name = name
    if is_mc:
        sample_name = f"{name}__sample_ctau_{sample_ctau:g}__reweight_ctau_{reweight_ctau:g}"
    sample_dir = outdir / cutflow_tag / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)

    cutflow_df.to_csv(sample_dir / "cutflow.csv", index=False)
    abcd_df.to_csv(sample_dir / "abcd.csv", index=False)
    cutflow_df.to_string(sample_dir / "cutflow.txt", index=False)
    abcd_df.to_string(sample_dir / "abcd.txt", index=False)

    

def main() -> None:
    args = parse_args()

    sample_defs = build_sample_defs(args)
    

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for name, cfg in sample_defs.items():
        run_sample(
            name,
            cfg,
            args.cutflow,
            args.size_cut,
            args.dphi_cut,
            args.blind,
            outdir,
            sample_ctau_override=args.sample_ctau,
            reweight_ctau_override=args.reweight_ctau,
        )


if __name__ == "__main__":
    main()
