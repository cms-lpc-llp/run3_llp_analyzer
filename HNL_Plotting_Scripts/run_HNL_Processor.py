#!/usr/bin/env python3

import argparse
import os
import sys
from typing import List
import warnings
import numpy as np

helper_dir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "python", "HNL_Plotting_HelperFunctions")
)
sys.path.insert(0, helper_dir)

warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    module="coffea\\.nanoevents\\.schemas\\..*",
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
warnings.filterwarnings(
    "ignore",
    message="Port 8787 is already in use.*",
    category=UserWarning,
    module="distributed\\.node",
)
warnings.filterwarnings(
    "ignore",
    message="Creating scratch directories is taking a surprisingly long time.*",
    category=UserWarning,
    module="contextlib",
)

import MuonSystemReader
import analysis_helpers
import HNL_Processor_v2
import HNL_Processor_v2_e
import HNL_Processor_v2_mu
import Processing_Helpers
import matplotlib.pyplot as plt


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
        "--not_normalized",
        action="store_true",
        help="Option to not normalize the histograms plotted (default: off for speed).",
    )
    parser.add_argument(
        "--from-tau",
        action="store_true",
        help="Apply fromTau option for e/mu processors (if supported).",
    )
    parser.add_argument(
        "--cut-based-id",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Use cut-based electron ID (electron only).",
    )
    parser.add_argument(
        "--iso",
        action="append",
        default=None,
        help="Isolation WP to apply (e or muon only).",
    )
    parser.add_argument(
        "--linscale",
        action="store_true",
        help="Plot Everything in Linear Scale",
    )
    parser.add_argument(
        "--fine-only",
        action="store_true",
        help="Only run fine-binning plots (data-only pass).",
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
    if args.flavor != "e" and args.cut_based_id is not None:
        raise ValueError(
            "--cut-based-id/--no-cut-based-id is only valid for --flavor e."
        )
    if args.flavor != "mu" and args.flavor != "e" and args.iso is not None:
        raise ValueError("--iso is only valid for --flavor mu or e.")
    linscale = args.linscale
    print(linscale)
    hists_to_process = None
    if args.hists:
        hists_to_process = [h.strip() for h in args.hists.split(",") if h.strip()]

    is_mc_list = to_bool_list(args.is_mc)

    iso_list = None
    if args.flavor == "mu" or args.flavor == "e":
        if args.iso is None:
            iso_list = ["tight"] * len(args.sample)
        else:
            if len(args.iso) != len(args.sample):
                raise ValueError(
                    f"--iso must be provided once per sample "
                    f"(got {len(args.iso)}, need {len(args.sample)})."
                )
            iso_list = args.iso

    outputs = []
    data_outputs_fine = []
    data_labels = []
    fine_bin_multiplier = 10
    for i, (sample_path, label, is_mc) in enumerate(
        zip(args.sample, args.label, is_mc_list)
    ):
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
        if not args.fine_only:
            if args.flavor in ("e", "mu"):
                process_kwargs = dict(
                    events=events,
                    hists_to_process=hists_to_process_local,
                    fromTau=args.from_tau,
                )
                if args.flavor == "e":
                    process_kwargs["cut_based_ID"] = (
                        True if args.cut_based_id is None else args.cut_based_id
                    )
                if args.flavor == "mu" or args.flavor == "e":
                    process_kwargs["iso"] = iso_list[i]
                    print(
                        f"[plot-debug] processing {label} with iso={iso_list[i]}",
                        flush=True,
                    )
                output = processor.process(**process_kwargs)
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
                process_kwargs = dict(
                    events=events,
                    hists_to_process=hists_to_process,
                    fromTau=args.from_tau,
                )
                if args.flavor == "e":
                    process_kwargs["cut_based_ID"] = (
                        True if args.cut_based_id is None else args.cut_based_id
                    )
                if args.flavor == "mu" or args.flavor == "e":
                    process_kwargs["iso"] = iso_list[i]
                    print(
                        f"[plot-debug] fine processing {label} with iso={iso_list[i]}",
                        flush=True,
                    )
                output_fine = processor_fine.process(**process_kwargs)
            else:
                output_fine = processor_fine.process(
                    events,
                    hists_to_process=hists_to_process,
                )
            data_outputs_fine.append(output_fine)
            data_labels.append(label)

    fraction = (not args.not_normalized)
    if not args.fine_only:
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

        root_dicts = [
            analysis_helpers.hist_dict_to_fraction_dict(o, fraction=fraction)
            for o in outputs
        ]

        all_keys = set()
        for d in root_dicts:
            all_keys.update(d.keys())

        line = f"[plot-debug] total unique hist keys across samples: {len(all_keys)}"
        print(line, flush=True)
        debug_lines.append(line)
        debug_lines.append(f"[plot-debug] keys: {', '.join(sorted(all_keys))}")
        
        with open(os.path.join(outdir, "plot_debug_summary.txt"), "w") as debug_file:
            debug_file.write("\n".join(debug_lines))
        
        print(fraction)
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
                    logscale=(not linscale),
                    fraction=fraction,
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
                    logscale=(not linscale),
                    fraction=fraction,
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
                    logscale=(not linscale),
                    fraction=fraction,
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
                    logscale=(not linscale),
                    fraction=fraction,
                    outdir=outdir,
                    filename=key,
                )
            else:
                raise ValueError("Only up to 4 samples are supported for plotting.")

    if data_outputs_fine:
        root_dicts_fine = [
            analysis_helpers.hist_dict_to_fraction_dict(o, fraction=fraction)
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
                logscale=(not linscale),
                fraction=fraction,
                outdir=outdir_fine,
            )
        elif len(root_dicts_fine) == 2:
            analysis_helpers.plot_hists_2_dicts(
                root_dicts_fine[0],
                root_dicts_fine[1],
                data_labels[0],
                data_labels[1],
                logscale=(not linscale),
                fraction=fraction,
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
                logscale=(not linscale),
                fraction=fraction,
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
                logscale=(not linscale),
                fraction=fraction,
                outdir=outdir_fine,
            )
        else:
            raise ValueError("Only up to 4 data samples are supported for fine-binning plotting.")

    # Plot event counts per ISO working point (uses electron_pt/muon_pt counts).
    def plot_iso_event_counts(outputs_list, outdir, labels, flavor):
        if not outputs_list:
            return
        hist_key = "electron_pt" if flavor == "e" else "muon_pt"
        counts = []
        used_labels = []
        for label, out in zip(labels, outputs_list):
            if not isinstance(out, dict) or hist_key not in out:
                continue
            hist_obj = out[hist_key]
            count = float(np.sum(hist_obj.view()))
            counts.append(count)
            used_labels.append(label)
        if not counts:
            return
        os.makedirs(outdir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7, 4))
        x = np.arange(len(counts))
        ax.bar(x, counts, color="#2b6cb0")
        ax.set_xticks(x)
        ax.set_xticklabels(used_labels, rotation=30, ha="right")
        ax.set_ylabel("Event Count")
        ax.set_title(f"{flavor} iso working point event counts")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "iso_event_counts.png"))
        with open(os.path.join(outdir, "iso_event_counts.txt"), "w") as f:
            for label, count in zip(used_labels, counts):
                f.write(f"{label}\t{count}\n")
        plt.close(fig)

    if args.flavor in ("e", "mu"):
        if not args.fine_only:
            plot_iso_event_counts(outputs, outdir, args.label, args.flavor)
        elif data_outputs_fine:
            plot_iso_event_counts(data_outputs_fine, outdir_fine, data_labels, args.flavor)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
