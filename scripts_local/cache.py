from __future__ import annotations
from multiprocessing import Pool
from pathlib import Path

import click
import numpy as np
import uproot as up
import h5py as h5
from tqdm import tqdm
import os

def load(fname: str | Path) -> np.ndarray:
    try:
        rm_flag = False
        if "cmsxrootd.fnal.gov" in fname:
            if Path(fname.replace("root://cmsxrootd.fnal.gov//","/storage/cms/")).exists():
                fname_temp = fname.replace("root://cmsxrootd.fnal.gov//","/storage/cms/")
            else:
                os.system(f"xrdcp {fname} .")
                fname_temp = str(Path(fname).name)
                rm_flag = True
        else:
            fname_temp = fname
        with up.open(fname_temp) as f:  # type: ignore
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
        if rm_flag: os.system(f"rm {fname_temp}")
        return np.asarray(np.stack([run, lumi, event], axis=1))
    except Exception as e:
        print(f'Error: {e} in {fname}')
        return np.array([])


def save_cache(fname: str | Path, data: dict[str, np.ndarray]) -> None:
    with h5.File(fname, 'w') as f:
        for k, v in tqdm(data.items(), desc='Saving'):
            f.create_dataset(k, data=v, compression='lzf')


def to_data(idxs: list[np.ndarray], paths: list[Path]) -> dict[str, np.ndarray]:
    data = {}
    for arr, p in zip(idxs, paths):
        name = p.with_suffix('').relative_to(p.parents[5]).__str__().replace('/', '=')
        data[name] = arr
    return data


def batch_load(base_path: Path | str, j: int | None = None):
    base_path = Path(base_path)
    pool = Pool(j)
    if base_path.is_dir():
        paths = list(base_path.glob('**/*.root'))
    else:
        paths = base_path.read_text().splitlines()
    r = pool.imap(load, paths)
    idxs = list(tqdm(r, total=len(paths), desc='Loading'))
    paths = [Path(p) for p in paths]
    return to_data(idxs, paths)


@click.command()
@click.option('--inp', '-i', type=str, required=True)
@click.option('--out', '-o', type=str, required=True)
@click.option('--jobs', '-j', type=int, default=None)
def main(inp: str, out: str, jobs: int | None):

    print(inp, out, jobs)
    data = batch_load(inp, j=jobs)
    save_cache(out, data)


if __name__ == '__main__':
    main()
