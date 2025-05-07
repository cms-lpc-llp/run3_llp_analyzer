import uproot as up
from multiprocessing import Pool
import numpy as np
import awkward as ak
from pathlib import Path
import json
import argparse
from tqdm import tqdm


def load_run_lumi(path: Path):
    try:
        with up.open(path) as f:  # type: ignore
            run = f['Events/run'].array().to_numpy()  # type: ignore
            lumi = f['Events/luminosityBlock'].array().to_numpy()  # type: ignore
    except Exception as e:
        print(f"Error: {path} {e}")
        run, lumi = [], []

    run_lumi = np.stack([run, lumi], axis=1)
    return np.unique(run_lumi, axis=0)


def _run_lumi_to_dict(run_lumi: np.ndarray):
    ret = {}
    cur_run = -1
    lumi0, lumi1 = -1, -1

    for run, lumi in run_lumi:
        if run != cur_run or lumi != lumi1 + 1:
            ret.setdefault(cur_run, []).append((lumi0, lumi1))
            lumi0 = int(lumi)
        cur_run = int(run)
        lumi1 = int(lumi)
    ret.setdefault(str(cur_run), []).append((lumi0, lumi1))
    ret.pop(-1)
    return ret


def load_dict(root: str, nproc: int = 8):
    pools = Pool(nproc)
    if Path(root).is_dir():
        paths = list(Path(root).rglob("**/*.root"))
    else:
        with open(root) as f:
            paths = [Path(line.strip()) for line in f.readlines()]
    run_lumis = []
    r = pools.imap_unordered(load_run_lumi, paths)
    for run_lumi in tqdm(r, total=len(paths)):
        run_lumis.append(run_lumi)
    pools.close()
    pools.join()
    run_lumis_ = np.unique(np.concatenate(run_lumis, axis=0), axis=0)
    return _run_lumi_to_dict(run_lumis_)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, required=True)
    parser.add_argument('--output', '-o', type=str, required=True)
    parser.add_argument('--jobs', '-j', type=int, default=8)
    args = parser.parse_args()

    run_lumi_dict = load_dict(args.input, args.jobs)
    with open(args.output, 'w') as f:
        json.dump(run_lumi_dict, f)
