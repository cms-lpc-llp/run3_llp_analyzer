import h5py as h5
import numpy as np
from pathlib import Path
import uproot as up
import argparse
from glob import glob
from multiprocessing import Pool
from tqdm import tqdm

from find_match import main as find_match

def batch(l:list, n:int):
    return [l[i:i+n] for i in range(0, len(l), n)]

def worker(args):
    n = find_match(*args[:-1])
    return n, args[-1]

def count(cache:str|Path):
    with h5.File(cache, 'r') as f:
        return len(f['cache'])

def main(path: str, batch_ntuple: int, jobs: int, output: str):
    cache_files = glob(f'{path}/**/*.h5', recursive=True)
    nano_files_by_era: dict[str, list[str]] = {}
    ntuple_files_by_era: dict[str, list[str]] = {}
    for fpath in cache_files:
        i=fpath.index('Run')
        era = fpath[i:i+8]
        fname = Path(fpath).name
        if fname.startswith('NANO'):
            nano_files_by_era.setdefault(era, []).append(fpath)
        else:
            assert 'ntupler' in fname
            ntuple_files_by_era.setdefault(era, []).append(fpath)
    
    for era in nano_files_by_era.keys():
        print(f'Era: {era} - {len(nano_files_by_era[era]):<5} nanos, {len(ntuple_files_by_era[era]):<5} ntuples')
    print(f'Era:  Total   - {sum([len(nano_files_by_era[era]) for era in nano_files_by_era.keys()]):<5} nanos, {sum([len(ntuple_files_by_era[era]) for era in ntuple_files_by_era.keys()]):<5} ntuples')

    args = []
    for era in sorted(nano_files_by_era.keys()):
        nanos = nano_files_by_era[era]
        ntuples = ntuple_files_by_era[era]
        for i, _ntuples in enumerate(batch(ntuples, batch_ntuple)):
            for nano in nanos:
                _out = Path(output) / era / f'{Path(nano).name}_{i}.h5'
                args.append(([nano], _ntuples, _out, era))
    
    n_matchs_per_era:dict[str,int] = {}

    pool = Pool(jobs)
    with tqdm(total=len(args)) as pbar:
        for r in pool.imap_unordered(worker, args):
            n_matchs_per_era[r[1]] = n_matchs_per_era.get(r[1], 0) + r[0]
            pbar.update(1)
    
    for era in n_matchs_per_era.keys():
        print(f'Era: {era} - {n_matchs_per_era[era]} matches')

    n_available_per_era:dict[str,int] = {}
    for k,v in ntuple_files_by_era.items():
        n_available_per_era[k] = sum([count(f) for f in v])

    for era in sorted(n_available_per_era.keys()):
        nm = n_matchs_per_era.get(era, 0)
        na = n_available_per_era[era]
        print(f'Era: {era} - {nm:<8} / {na:<8} ({1-nm/(na+1e-9):.2%} lost)')
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str)
    parser.add_argument('-b', '--batch-ntuple', type=int, required=False, default=999999)
    parser.add_argument('-j', '--jobs', type=int, required=False, default=1)
    parser.add_argument('-o', '--output', type=str, required=True)
    args = parser.parse_args()

    main(args.path, args.batch_ntuple, args.jobs, args.output)


