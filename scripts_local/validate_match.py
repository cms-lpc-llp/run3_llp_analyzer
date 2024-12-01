import h5py as h5
import numpy as np
from pathlib import Path
import uproot as up
import argparse


def load_rle(fs: list[str | Path]):
    rles = []
    for f in fs:
        with up.open(f) as f:
            if 'ntuples;1' in f.keys():
                runs = f['ntuples/llp/runNum'].array().to_numpy()
                events = f['ntuples/llp/eventNum'].array().to_numpy()
                lumis = f['ntuples/llp/lumiNum'].array().to_numpy()
            else:
                runs = f['Events/run'].array().to_numpy()
                events = f['Events/event'].array().to_numpy()
                lumis = f['Events/luminosityBlock'].array().to_numpy()
        _rle = np.stack([runs, lumis, events], axis=1)
        rles.append(_rle)
    return np.concatenate(rles)


def main(cache_file: str | Path):
    with up.open(cache_file) as f:
        nano_files = [str(x) for x in f['inp1/paths'].array()]
        ntuple_files = [str(x) for x in f['inp2/paths'].array()]
        l1s = f['inp1/lengths'].array().to_numpy()
        l2s = f['inp2/lengths'].array().to_numpy()
        idxs_nano = f['idx/1'].array().to_numpy()
        idxs_ntuple = f['idx/2'].array().to_numpy()
        rle1 = load_rle(nano_files)
        rle2 = load_rle(ntuple_files)
        # print(len(rle1), len(rle2))
        # print(l1s, l2s)
        # for r, l, e in rle2:
        #     print(f'{r}:{l}:{e}')

        for i1, i2 in zip(idxs_nano, idxs_ntuple):
            print(f'{i1} {i2}')
            # print(rle1[i1], '-', rle2[i2])
            assert np.all(rle1[i1] == rle2[i2])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)
    args = parser.parse_args()

    main(args.file)
