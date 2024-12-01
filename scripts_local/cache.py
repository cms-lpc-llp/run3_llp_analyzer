from multiprocessing import Pool
from pathlib import Path

import click
import numpy as np
import uproot as up
import h5py as h5
from tqdm import tqdm


def load(fname: str | Path) -> np.ndarray:
    try:
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
    except Exception as e:
        print(f'Error: {e} in {fname}')
        return np.array([])


def save_cache(fname: str | Path, data: dict[str, np.ndarray]) -> None:
    with h5.File(fname, 'w') as f:
        for k, v in tqdm(data.items(), desc='Saving'):
            f.create_dataset(k, data=v, compression='lzf')


def to_data(idxs: list[np.ndarray], paths: list[str]) -> dict[str, np.ndarray]:
    data = {}
    for arr, p in zip(idxs, paths):
        data[p] = arr
    return data


def batch_load(base_path: Path | str, j: int | None = None):
    base_path = Path(base_path)
    pool = Pool(j)
    if base_path.is_dir():
        paths = [str(p) for p in base_path.glob('**/*.root')]
    else:
        paths = base_path.read_text().splitlines()
    r = pool.imap(load, paths)
    idxs = list(tqdm(r, total=len(paths), desc='Loading'))
    return to_data(idxs, paths)


@click.command()
@click.option('--inp', '-i', type=str, required=True)
@click.option('--out', '-o', type=str, required=True)
@click.option('--jobs', '-j', type=int, default=None)
def main(inp: str, out: str, jobs: int | None):
    data = batch_load(inp, j=jobs)
    save_cache(out, data)


if __name__ == '__main__':
    main()
