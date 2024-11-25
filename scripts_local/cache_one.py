from multiprocessing import Pool
from pathlib import Path

import argparse
import numpy as np
import uproot as up
import h5py as h5
from tqdm import tqdm


def load(fname: str | Path) -> np.ndarray:
    with up.open(fname) as f:  # type: ignore
        keys = f.keys()
        if 'ntuples;1' in keys:
            run = f['ntuples/llp/runNum'].array()  # type: ignore
            lumi = f['ntuples/llp/lumiNum'].array()  # type: ignore
            event = f['ntuples/llp/eventNum'].array()  # type: ignore
        elif 'Events;1' in keys:
            run = f['Events/run'].array()  # type: ignore
            lumi = f['Events/luminosityBlock'].array()  # type: ignore
            event = f['Events/event'].array()  # type: ignore
        else:
            raise ValueError(f'Could not find ntuple or Events in {fname}')
    return np.asarray(np.stack([run, lumi, event], axis=1))


def save_cache(fname: str | Path, key: str, data: np.ndarray) -> None:
    with h5.File(fname, 'w') as f:
        f.create_dataset('cache', data=data, compression='lzf')
        f.attrs['key'] = key


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inp', '-i', type=str, required=True)
    parser.add_argument('--out', '-o', type=str, required=True)
    args = parser.parse_args()

    inp = args.inp
    out = args.out
    data = load(inp)
    save_cache(out, inp, data)
    print(f'Saved {out}')
