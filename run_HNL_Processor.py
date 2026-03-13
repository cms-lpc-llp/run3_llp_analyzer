#!/usr/bin/env python3

import argparse
import os
import sys
from typing import List
import warnings

helper_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "python", "HNL_Plotting_HelperFunctions")
)
sys.path.insert(0, helper_dir)

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
import analysis_helpers
import HNL_Processor_v2
import HNL_Processor_v2_e
import HNL_Processor_v2_mu
import Processing_Helpers


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run HNL processors to produce ROOT-style plots (matching the "
            "Preliminary_Sensitivity_Plots notebooks)."
        )
    )
    parser.add_argument(
        "--flavor",
        required=True,
        choices=("e", "mu", "tau"),
        help="Lepton flavor to select the correct processor.",
    )
    parser.add_argument(
        "--sample",
        required=True,
        action="append",
        help="Input ROOT file path. Pass multiple times for multiple samples.",
    )
    parser.add_argument(
        "--label",
        required=True,
        action="append",
        help="Label for each sample. Must match the number of --sample args.",
    )
    parser.add_argument(
        "--is-mc",
        required=True,
        action="append",
        help="Whether each sample is MC (true/false). Must match --sample args.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory name (created under selection_variable_plots).",
    )
    parser.add_argument(
        "--hists",
        default=None,
        help="Comma-separated list of histogram names to process (optional).",
    )
    parser.add_argument(
        "--apply-gen-info",
        action="store_true",
        help="Enable gen-level histograms for MC (default: off for speed).",
    )
    parser.add_argument(
        "--from-tau",
        action="store_true",
        help="Apply fromTau option for e/mu processors (if supported).",
    )
    parser.add_argument(
        "--cut-based-id",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use cut-based electron ID (electron flavor only).",
    )
    return parser.parse_args()


def to_bool_list(values: List[str]) -> List[bool]:
    truthy = {"1", "true", "yes", "y", "t"}
    falsy = {"0", "false", "no", "n", "f"}
    out = []
    for value in values:
        v = value.strip().lower()
        if v in truthy:
            out.append(True)
        elif v in falsy:
            out.append(False)
        else:
            raise ValueError(f"Invalid boolean value for --is-mc: {value}")
    return out


def main() -> int:
    print(f"[plot-debug] run_HNL_Processor path: {__file__}", flush=True)
    args = parse_args()

    if not (len(args.sample) == len(args.label) == len(args.is_mc)):
        raise ValueError(
            "Lengths of --sample, --label, and --is-mc must match."
        )

    if args.flavor == "tau":
        processor_cls = HNL_Processor_v2.HNL_Processor_v2
        trigger = "tau"
    elif args.flavor == "e":
        processor_cls = HNL_Processor_v2_e.HNL_Processor_v2_e
        trigger = "electron"
    else:
        processor_cls = HNL_Processor_v2_mu.HNL_Processor_v2_mu
        trigger = "muon"

    hists_to_process = None
    if args.hists:
        hists_to_process = [h.strip() for h in args.hists.split(",") if h.strip()]

    is_mc_list = to_bool_list(args.is_mc)

    outputs = []
    data_outputs_fine = []
    data_labels = []
    fine_bin_multiplier = 10
    for sample_path, label, is_mc in zip(args.sample, args.label, is_mc_list):
        if is_mc:
            events = MuonSystemReader.loadTree_nanoFactory(
                sample_path, isMC=True, trigger=trigger
            )
        else:
            data_events_list = Processing_Helpers.make_data_list(sample_path, args.flavor)
            events = MuonSystemReader.loadTree_nanoFactory(
                data_events_list, isMC=False, trigger=trigger
            )
        processor = processor_cls(isMC=is_mc, applyGenInfo=args.apply_gen_info)
        hists_to_process_local = hists_to_process
        if args.flavor in ("e", "mu"):
            output = processor.process(
                events,
                hists_to_process=hists_to_process_local,
                fromTau=args.from_tau,
                cut_based_ID=args.cut_based_id if args.flavor == "e" else True,
            )
        else:
            output = processor.process(
                events,
                hists_to_process=hists_to_process_local,
            )
        outputs.append(output)
        if not is_mc:
            processor_fine = processor_cls(
                isMC=False,
                applyGenInfo=args.apply_gen_info,
                bin_multiplier=fine_bin_multiplier,
            )
            if args.flavor in ("e", "mu"):
                output_fine = processor_fine.process(
                    events,
                    hists_to_process=hists_to_process,
                    fromTau=args.from_tau,
                    cut_based_ID=args.cut_based_id if args.flavor == "e" else True,
                )
            else:
                output_fine = processor_fine.process(
                    events,
                    hists_to_process=hists_to_process,
                )
            data_outputs_fine.append(output_fine)
            data_labels.append(label)

    outdir = os.path.join(
        os.path.dirname(__file__),
        "selection_variable_plots",
        args.outdir,
    )
    os.makedirs(outdir, exist_ok=True)
    debug_lines = []
    for label, out in zip(args.label, outputs):
        if isinstance(out, dict):
            keys = sorted(out.keys())
            preview = ", ".join(keys[:5])
            more = f" (+{len(keys) - 5} more)" if len(keys) > 5 else ""
            line = f"[plot-debug] {label}: {len(keys)} hist keys -> {preview}{more}"
            print(line, flush=True)
            debug_lines.append(line)
        else:
            line = f"[plot-debug] {label}: output is {type(out)} (expected dict)"
            print(line, flush=True)
            debug_lines.append(line)

    root_dicts = [analysis_helpers.hist_dict_to_fraction_dict(o) for o in outputs]

    all_keys = set()
    for d in root_dicts:
        all_keys.update(d.keys())

    line = f"[plot-debug] total unique hist keys across samples: {len(all_keys)}"
    print(line, flush=True)
    debug_lines.append(line)
    debug_lines.append(f"[plot-debug] keys: {', '.join(sorted(all_keys))}")

    with open(os.path.join(outdir, "plot_debug_summary.txt"), "w") as debug_file:
        debug_file.write("\n".join(debug_lines))
    for key in sorted(all_keys):
        present = [i for i, d in enumerate(root_dicts) if key in d]
        if not present:
            continue
        hists = [root_dicts[i][key][0] for i in present]
        title = root_dicts[present[0]][key][1]
        x_label = root_dicts[present[0]][key][2]
        labels = [args.label[i] for i in present]

        if len(hists) == 1:
            analysis_helpers.plot_one_root_hists(
                hists[0],
                title,
                x_label,
                labels[0],
                logscale=True,
                fraction=True,
                outdir=outdir,
                filename=key,
            )
        elif len(hists) == 2:
            analysis_helpers.plot_two_root_hists(
                hists[0],
                hists[1],
                title,
                x_label,
                labels[0],
                labels[1],
                logscale=True,
                fraction=True,
                outdir=outdir,
                filename=key,
            )
        elif len(hists) == 3:
            analysis_helpers.plot_three_root_hists(
                hists[0],
                hists[1],
                hists[2],
                title,
                x_label,
                labels[0],
                labels[1],
                labels[2],
                logscale=True,
                fraction=True,
                outdir=outdir,
                filename=key,
            )
        elif len(hists) == 4:
            analysis_helpers.plot_four_root_hists(
                hists[0],
                hists[1],
                hists[2],
                hists[3],
                title,
                x_label,
                labels[0],
                labels[1],
                labels[2],
                labels[3],
                logscale=True,
                fraction=True,
                outdir=outdir,
                filename=key,
            )
        else:
            raise ValueError("Only up to 4 samples are supported for plotting.")

    if data_outputs_fine:
        root_dicts_fine = [
            analysis_helpers.hist_dict_to_fraction_dict(o)
            for o in data_outputs_fine
        ]
        outdir_fine = os.path.join(
            os.path.dirname(__file__),
            "selection_variable_plots",
            f"{args.outdir}_fineBinning",
        )

        if len(root_dicts_fine) == 1:
            analysis_helpers.plot_hists_one_dicts(
                root_dicts_fine[0],
                data_labels[0],
                logscale=True,
                fraction=True,
                outdir=outdir_fine,
            )
        elif len(root_dicts_fine) == 2:
            analysis_helpers.plot_hists_2_dicts(
                root_dicts_fine[0],
                root_dicts_fine[1],
                data_labels[0],
                data_labels[1],
                logscale=True,
                fraction=True,
                outdir=outdir_fine,
            )
        elif len(root_dicts_fine) == 3:
            analysis_helpers.plot_hists_3_dicts(
                root_dicts_fine[0],
                root_dicts_fine[1],
                root_dicts_fine[2],
                data_labels[0],
                data_labels[1],
                data_labels[2],
                logscale=True,
                fraction=True,
                outdir=outdir_fine,
            )
        elif len(root_dicts_fine) == 4:
            analysis_helpers.plot_hists_4_dicts(
                root_dicts_fine[0],
                root_dicts_fine[1],
                root_dicts_fine[2],
                root_dicts_fine[3],
                data_labels[0],
                data_labels[1],
                data_labels[2],
                data_labels[3],
                logscale=True,
                fraction=True,
                outdir=outdir_fine,
            )
        else:
            raise ValueError("Only up to 4 data samples are supported for fine-binning plotting.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
