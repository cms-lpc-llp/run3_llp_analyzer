from multiprocessing import Pool
from pathlib import Path

import click
import numpy as np
import uproot as up
import h5py as h5
from tqdm import tqdm


def save_cache(fname: str | Path, data: np.ndarray):
    with up.recreate(fname) as f:
        tree = {'run': data[:, 0].astype(np.uint32), 'event': data[:, -1]}
        f['events'] = tree


def worker(args):
    save_cache(*args)


@click.command()
@click.option('--inp', '-i', type=str, required=True)
@click.option('--out', '-o', type=str, required=True)
@click.option('--jobs', '-j', type=int, default=None)
def main(inp: str, out: str, jobs: int | None):

    outp = Path(out)
    outp.mkdir(exist_ok=True, parents=True)
    with h5.File(inp, 'r') as f:
        l = len(f)

    def gen():
        with h5.File(inp, 'r') as f:
            for k, v in f.items():
                if v.shape == (0,):
                    continue
                yield outp / f'{k}.root', v[:]

    pool = Pool(jobs)
    r = pool.imap(worker, gen())
    list(tqdm(r, desc='Processing', total=l))


if __name__ == '__main__':
    main()
