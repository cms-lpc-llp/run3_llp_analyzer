#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Loop over signal reweight ctaus, run Run_Cutflow_ABCD for each point, "
            "then write datacards using data counts from a provided abcd.txt."
        )
    )
    parser.add_argument(
        "--cutflow",
        required=True,
        help="Cutflow YAML path relative to cuts_config, e.g. mu_SR/cuts_sensitivity_mu_DNN.yaml",
    )
    parser.add_argument(
        "--signal-name",
        action="append",
        required=True,
        help=(
            "Signal process name used in datacards (modelName/process label). "
            "Can be provided multiple times for multiple samples."
        ),
    )
    parser.add_argument(
        "--sample",
        action="append",
        required=True,
        help="Signal sample key/name passed to Run_Cutflow_ABCD. Can be repeated.",
    )
    parser.add_argument(
        "--sample-config",
        default=None,
        help="Optional sample config YAML passed to Run_Cutflow_ABCD.",
    )
    parser.add_argument(
        "--data-abcd",
        required=True,
        help="Path to data abcd.txt used to build list_data/background rates.",
    )
    parser.add_argument(
        "--ctaus",
        default="10,25,50,100,250,500,750,1000,2500,5000,7500,10000",
        help="Comma-separated reweight ctaus (mm).",
    )
    parser.add_argument(
        "--sample-ctau",
        action="append",
        type=float,
        default=None,
        help=(
            "Generated/sample ctau (mm) of the input signal sample. "
            "Provide once to apply to all samples, or once per --sample."
        ),
    )
    parser.add_argument(
        "--mass",
        type=float,
        required=True,
        help="Signal mass used in output datacard filename.",
    )
    parser.add_argument(
        "--outdir",
        default="SR_Cutflows",
        help="Output directory for Run_Cutflow_ABCD.",
    )
    parser.add_argument(
        "--datacards-outdir",
        required=True,
        help="Directory where datacards will be written.",
    )
    parser.add_argument(
        "--signal-region",
        default="A",
        help="Datacard signal region label (kept for compatibility).",
    )
    parser.add_argument(
        "--prefix",
        default="bin",
        help="Datacard ABCD nuisance prefix.",
    )
    parser.add_argument(
        "--sig-unc",
        type=float,
        default=0.2,
        help="Flat lnN signal uncertainty for each ABCD bin.",
    )
    parser.add_argument(
        "--sig-unc-name",
        default="dummy",
        help="Name of the flat signal uncertainty nuisance.",
    )
    parser.add_argument(
        "--size-cut",
        type=float,
        default=160.0,
        help="Pass-through to Run_Cutflow_ABCD.",
    )
    parser.add_argument(
        "--dphi-cut",
        type=float,
        default=2.0,
        help="Pass-through to Run_Cutflow_ABCD.",
    )
    parser.add_argument(
        "--blind",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Pass-through to Run_Cutflow_ABCD.",
    )
    return parser.parse_args()


def cutflow_tag_from_path(cutflow_cfg: str) -> str:
    tag = cutflow_cfg.strip().replace("\\", "/")
    if tag.endswith(".yaml"):
        tag = tag[:-5]
    tag = tag.replace("/", "__")
    return tag or "cutflow"


def parse_ctaus(value: str) -> List[float]:
    out = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        out.append(float(token))
    if not out:
        raise ValueError("No ctaus were provided.")
    return out


def _ctau_str(value: float) -> str:
    return f"{value:g}"


def resolve_sample_name_pairs(samples: List[str], signal_names: List[str]) -> List[tuple[str, str]]:
    if len(samples) == len(signal_names):
        return list(zip(samples, signal_names))
    raise ValueError(
        "Mismatched number of --sample and --signal-name arguments. "
        "Provide one --signal-name per --sample."
    )


def resolve_sample_ctaus(samples: List[str], sample_ctaus: Optional[List[float]]) -> List[float]:
    if not sample_ctaus:
        return [1000.0] * len(samples)
    if len(sample_ctaus) == 1:
        return [sample_ctaus[0]] * len(samples)
    if len(sample_ctaus) == len(samples):
        return sample_ctaus
    raise ValueError(
        "Mismatched number of --sample-ctau arguments. "
        "Provide one value for all samples or one per --sample."
    )


def parse_abcd_txt(path: Path) -> Dict[str, float]:
    lines = path.read_text().splitlines()
    counts_by_bin: Dict[str, float] = {}
    for line in lines:
        m_bin = re.search(r"bin\s+([ABCD])", line, flags=re.IGNORECASE)
        if not m_bin:
            continue
        m_count = re.search(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*\+\-", line)
        if not m_count:
            continue
        counts_by_bin[m_bin.group(1).upper()] = float(m_count.group(1))
    missing = [b for b in ["A", "B", "C", "D"] if b not in counts_by_bin]
    if missing:
        raise ValueError(f"Missing bins {missing} while parsing {path}")
    return counts_by_bin


def to_datacard_rate_order(counts_by_bin: Dict[str, float]) -> List[float]:
    # Match historical notebook convention: [D, B, A, C].
    return [counts_by_bin["D"], counts_by_bin["B"], counts_by_bin["A"], counts_by_bin["C"]]


def make_datacard_2tag(
    out_datacards_dir: Path,
    model_name: str,
    signal_rate: Dict[str, List[float]],
    norm: float,
    bkg_rate: List[float],
    observation: List[float],
    sig_unc: Dict[str, List[List[float]]],
    sig_unc_name: List[str],
    signal_region: str,
    prefix: str,
    mass: float,
    ctau: float,
) -> Path:
    a, b, c, d = bkg_rate[0], bkg_rate[1], bkg_rate[2], bkg_rate[3]
    n_sig = len(signal_rate.keys())
    model_dir = out_datacards_dir / model_name
    model_dir.mkdir(parents=True, exist_ok=True)
    card_path = model_dir / f"mN{mass:g}_ctau{ctau:g}_norm{norm:.5f}.txt"

    with card_path.open("w") as text_file:
        text_file.write(f"# signal norm {norm} \n")
        text_file.write("imax 4 \n")
        text_file.write(f"jmax {n_sig} \n")
        text_file.write("kmax * \n")
        text_file.write("shapes * * FAKE \n")
        text_file.write("--------------- \n")
        text_file.write("--------------- \n")
        text_file.write("bin \t chA \t chB \t chC \t chD \n")
        text_file.write(
            "observation \t {0:6.2f} \t {1:6.2f} \t {2:6.2f} \t {3:6.2f} \n".format(
                observation[0], observation[1], observation[2], observation[3]
            )
        )
        text_file.write("------------------------------ \n")
        text_file.write("bin " + "\t chA " * (1 + n_sig) + "\t chB " * (1 + n_sig) + "\t chC " * (1 + n_sig) + "\t chD " * (1 + n_sig) + "\n")
        process_name = "\t " + (" \t ").join(list(signal_rate.keys())) + "\t bkg "
        text_file.write("process " + process_name * 4 + "\n")
        process_number = "\t " + (" \t ").join(list((np.arange(n_sig) * -1).astype(str))) + "\t 1"
        text_file.write("process " + process_number * 4 + "\n")

        rate_string = "rate"
        for i in range(4):
            for _, v in signal_rate.items():
                rate_string += "\t {0:e} ".format(v[i])
            rate_string += "\t 1 "
        text_file.write(rate_string + "\n")
        text_file.write("------------------------------ \n")

        text_file.write(prefix + "A   rateParam       chA     bkg     (@0*@2/@1)     " + prefix + "B," + prefix + "C," + prefix + "D \n")
        text_file.write(prefix + "B   rateParam       chB     bkg     {0:.2f}        [0,{1:.2f}] \n".format(b, c * 7 if b == 0 else b * 7))
        text_file.write(prefix + "C   rateParam       chC     bkg     {0:.2f}        [0,{1:.2f}] \n".format(c, c * 7))
        text_file.write(prefix + "D   rateParam       chD     bkg     {0:.2f}        [0,{1:.2f}] \n".format(d, c * 7 if d == 0 else d * 7))

        for k in signal_rate.keys():
            text_file.write("norm rateParam * {0} {1:.2f}  \n".format(k, norm))

        for k, v in sig_unc.items():
            assert len(sig_unc_name) == len(v)

        for i in range(len(sig_unc_name)):
            unc_text = sig_unc_name[i] + " \t lnN"
            for j in range(4):
                for _, v in sig_unc.items():
                    if v[i][j] == 0.0:
                        unc_text += " \t -"
                    else:
                        unc_text += " \t " + str(v[i][j] + 1)
                unc_text += "\t - "
            text_file.write(unc_text + " \n")

        # Keep arg for compatibility with existing workflows.
        _ = signal_region

    return card_path


def run_cutflow_for_ctau(
    run_script: Path,
    cutflow: str,
    sample: str,
    sample_config: str,
    sample_ctau: float,
    reweight_ctau: float,
    size_cut: float,
    dphi_cut: float,
    blind: bool,
    outdir: Path,
) -> None:
    cmd = [
        sys.executable,
        str(run_script),
        "--cutflow",
        cutflow,
        "--sample",
        sample,
        "--sample-ctau",
        _ctau_str(sample_ctau),
        "--reweight-ctau",
        _ctau_str(reweight_ctau),
        "--size-cut",
        _ctau_str(size_cut),
        "--dphi-cut",
        _ctau_str(dphi_cut),
        "--outdir",
        str(outdir),
    ]
    if sample_config:
        cmd.extend(["--sample-config", sample_config])
    if blind:
        cmd.append("--blind")
    else:
        cmd.append("--no-blind")
    subprocess.run(cmd, check=True)


def find_signal_abcd_path(outdir: Path, cutflow: str, sample: str, sample_ctau: float, reweight_ctau: float) -> Path:
    tag = cutflow_tag_from_path(cutflow)
    expected_dir = outdir / tag / f"{sample}__sample_ctau_{sample_ctau:g}__reweight_ctau_{reweight_ctau:g}"
    expected_file = expected_dir / "abcd.txt"
    if expected_file.exists():
        return expected_file

    fallback_pattern = f"{sample}__sample_ctau_{sample_ctau:g}__reweight_ctau_{reweight_ctau:g}"
    candidates = list((outdir / tag).glob(f"{fallback_pattern}/abcd.txt"))
    if len(candidates) == 1:
        return candidates[0]
    if not candidates:
        raise FileNotFoundError(f"Could not find signal abcd.txt for ctau={reweight_ctau:g} in {(outdir / tag)}")
    raise RuntimeError(f"Multiple matching signal abcd.txt files found for ctau={reweight_ctau:g}: {candidates}")


def main() -> None:
    args = parse_args()
    run_script = Path(__file__).resolve().parent / "Run_Cutflow_ABCD.py"
    if not run_script.exists():
        raise FileNotFoundError(f"Could not find {run_script}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data_abcd_path = Path(args.data_abcd)
    if not data_abcd_path.exists():
        raise FileNotFoundError(f"Data abcd file not found: {data_abcd_path}")

    ctaus = parse_ctaus(args.ctaus)
    data_counts = parse_abcd_txt(data_abcd_path)
    list_data = to_datacard_rate_order(data_counts)

    sig_unc = {"dummy": [[args.sig_unc, args.sig_unc, args.sig_unc, args.sig_unc]]}
    sig_unc_name = [args.sig_unc_name]

    out_datacards_dir = Path(args.datacards_outdir)
    out_datacards_dir.mkdir(parents=True, exist_ok=True)

    sample_signal_pairs = resolve_sample_name_pairs(args.sample, args.signal_name)
    sample_ctaus = resolve_sample_ctaus(args.sample, args.sample_ctau)

    for (sample, signal_name), sample_ctau in zip(sample_signal_pairs, sample_ctaus):
        print(f"Processing sample={sample} (sample_ctau={sample_ctau:g}) -> signal_name={signal_name}")
        for ctau in ctaus:
            print(f"Running ctau={ctau:g}")
            run_cutflow_for_ctau(
                run_script=run_script,
                cutflow=args.cutflow,
                sample=sample,
                sample_config=args.sample_config,
                sample_ctau=sample_ctau,
                reweight_ctau=ctau,
                size_cut=args.size_cut,
                dphi_cut=args.dphi_cut,
                blind=args.blind,
                outdir=outdir,
            )

            signal_abcd_path = find_signal_abcd_path(
                outdir=outdir,
                cutflow=args.cutflow,
                sample=sample,
                sample_ctau=sample_ctau,
                reweight_ctau=ctau,
            )
            signal_counts = parse_abcd_txt(signal_abcd_path)
            list_signal = to_datacard_rate_order(signal_counts)

            signal_rate = {signal_name: list_signal}
            if signal_rate[signal_name][0] == 0:
                raise ZeroDivisionError(
                    f"Signal rate in datacard A-bin is zero for sample={sample}, ctau={ctau:g}; cannot compute norm."
                )
            norm = 1.0 / signal_rate[signal_name][0]
            bkg_rate = list_data
            observation = bkg_rate

            card_path = make_datacard_2tag(
                out_datacards_dir=out_datacards_dir,
                model_name=signal_name,
                signal_rate=signal_rate,
                norm=norm,
                bkg_rate=bkg_rate,
                observation=observation,
                sig_unc=sig_unc,
                sig_unc_name=sig_unc_name,
                signal_region=args.signal_region,
                prefix=args.prefix,
                mass=args.mass,
                ctau=ctau,
            )
            print(f"Wrote datacard: {card_path}")


if __name__ == "__main__":
    main()
