#!/usr/bin/env python3

import os
from functools import lru_cache

import pandas as pd


LOOKUP_DIR = os.path.join(
    os.environ["CMSSW_BASE"],
    "src",
    "run3_llp_analyzer",
    "python",
    "HNL_Plotting_HelperFunctions",
    "modeling_lookup_tables",
)

LOOKUP_FILES = {
    "cluster_size": "Cluster_Size_Data_MC_Comparison_072426.csv",
    "dnn": "DNN_Data_MC_Comparison_072426.csv",
}


@lru_cache(maxsize=None)
def _load_lookup_table(cut_kind):
    lookup_path = os.path.join(LOOKUP_DIR, LOOKUP_FILES[cut_kind])
    table = pd.read_csv(lookup_path, usecols=["Data Cut", "MC Cut"])
    return table.astype(float)


def _lookup_mc_cut(data_cut, cut_kind):
    table = _load_lookup_table(cut_kind)
    idx = (table["Data Cut"] - float(data_cut)).abs().idxmin()
    matched_data_cut = float(table.loc[idx, "Data Cut"])
    mc_cut = float(table.loc[idx, "MC Cut"])
    return matched_data_cut, mc_cut


def _should_map_cut_value(value, cut_kind):
    if value is None:
        return False
    numeric_value = float(value)
    if cut_kind == "dnn":
        return 0.0 <= numeric_value < 1.0
    if cut_kind == "cluster_size":
        return 0.0 < numeric_value < 1000.0
    return False


def _format_cut_value(value):
    numeric_value = float(value)
    if numeric_value.is_integer():
        return str(int(numeric_value))
    return f"{numeric_value:.6g}"


def remap_cut_info_for_mc(cut_key, cut_info):
    branch = cut_info.get("branch")
    if branch == "cscRechitClusterDNN_bkgMC_plusBeamHalo":
        cut_kind = "dnn"
    elif branch in {"cscRechitClusterSize", "cscRechitCluster3Size"}:
        cut_kind = "cluster_size"
    else:
        return cut_info

    remapped_cut_info = dict(cut_info)
    for bound_key in ("minVal", "maxVal"):
        bound_value = remapped_cut_info.get(bound_key)
        if not _should_map_cut_value(bound_value, cut_kind):
            continue
        matched_data_cut, mc_cut = _lookup_mc_cut(bound_value, cut_kind)
        remapped_cut_info[bound_key] = mc_cut
        print(
            f"[MC cut lookup] {cut_key} ({branch}) {bound_key}: "
            f"data {_format_cut_value(bound_value)} -> MC {_format_cut_value(mc_cut)} "
            f"(matched data {_format_cut_value(matched_data_cut)})"
        )

    return remapped_cut_info


def remap_abcd_size_cut_for_mc(size_cut, *, context):
    matched_data_cut, mc_cut = _lookup_mc_cut(size_cut, "cluster_size")
    print(
        f"[MC cut lookup] {context} sizeCut: "
        f"data {_format_cut_value(size_cut)} -> MC {_format_cut_value(mc_cut)} "
        f"(matched data {_format_cut_value(matched_data_cut)})"
    )
    return mc_cut
