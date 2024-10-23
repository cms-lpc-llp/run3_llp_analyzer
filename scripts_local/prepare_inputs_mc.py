from __future__ import annotations
from pathlib import Path
import uproot as up
from multiprocessing import Pool
import h5py as h5
from tqdm import tqdm
import numpy as np
import json


def load_key(k: str, fname: str|Path) -> np.ndarray:
    with h5.File(fname, 'r') as f:
        arr = f[k][()]
        return arr

def worker(x:tuple[str, str|Path]):
    key, fname = x
    return key, load_key(key, fname)

def era_from_path(path:Path) -> str:
    if 'NanoAOD' in str(path):
        return "MC_"+str(path)[str(path).find("Summer"):str(path).find("NanoAOD")]
    else:
        ll = str(path.parents[2]).split('/')
        for l in ll: 
            if 'MC' in l:return l
        return ''
            # return path.parents[2].name.split('_')[5][3:8]
            # return str(path.parents[2])[index:index]


base_path = '/storage/af/user/christiw/login-1/christiw/LLP/Run3/CMSSW_10_6_30'
ntuple_base_path = Path(base_path + '/src/run3_llp_analyzer/lists/displacedJetMuonNtuple/V1p19/')
nano_base_path = Path(base_path + '/src/run3_llp_analyzer/lists/nanoAOD/')



campaigns = { # keyword for ntuples and nanoAOD
#   'MC_Summer22': 'ggH',
#   'MC_Summer22EE': 'ggH',
    'MC_Summer23'  : 'ggH',
   'MC_Summer23BPix': 'ggH',
    
}

for key, prod in campaigns.items():
    print(f"RUNNING FOR {key}")
    
    ntuple_list_path_all = list(ntuple_base_path.glob(f"*{key}/v2/*/{prod}*"))
    # nano_list_path = list(nano_base_path.glob(f"**/{k}/{v}*"))
    for ntuple_list_path in ntuple_list_path_all:

        nano_list_path = list(nano_base_path.glob(f"**/{key}/{ntuple_list_path.stem}*"))
        
        print(ntuple_list_path)
        print(nano_list_path)
        assert(len(nano_list_path)==1)
        # exit()
        # ntuple_list_path = ntuple_list_path[0]
        nano_list_path = nano_list_path[0]

        ntuple_cache_path = str(ntuple_list_path).replace('lists','/data/cache/').split('.')[0]
        nano_cache_path = str(nano_list_path).replace('lists','/data/cache/').split('.')[0]

        ntuple_files = ntuple_list_path.read_text().splitlines()
        nano_files = nano_list_path.read_text().splitlines()

        ntuple_files = [Path(p) for p in ntuple_files]
        nano_files = [Path(p) for p in nano_files]
        ntuple_k_map = {str(p.with_suffix('').relative_to(p.parents[5]).__str__().replace('/', '=')):p for p in ntuple_files}
        nano_k_map = {str(p.with_suffix('').relative_to(p.parents[5]).__str__().replace('/', '=')):p for p in nano_files}

        with h5.File(ntuple_cache_path, 'r') as f:
            ntuple_keys = list(f.keys())

        with h5.File(nano_cache_path, 'r') as f:
            nano_keys = list(f.keys())

        pool = Pool(4)
        r = pool.imap(worker, [(k, ntuple_cache_path) for k in ntuple_keys])
        ntuple = {}
        for k,arr in tqdm(r, total=len(ntuple_keys), desc='Loading NTuples'):
            path = ntuple_k_map[k]
            era = era_from_path(path)
            ntuple.setdefault(era, {})[k] = set(tuple(x) for x in arr)
        r = pool.imap(worker, [(k, nano_cache_path) for k in nano_keys])
        nano = {}
        for k,arr in tqdm(r, total=len(nano_keys), desc='Loading Nano'):
            # path = nano_k_map[k.split('_')[-1]]
            # print(path)
            path = nano_k_map[k]
            era = era_from_path(path)
            nano.setdefault(era, {})[k] = arr
        n_events = {k:[sum(len(vv) for vv in v.values()),0] for k,v in ntuple.items()}
        matchs = {}
        
        with tqdm(total=sum(len(x) for x in nano.values())) as pbar:
            for era, d_nano in nano.items():
                if era not in ntuple:continue
                d_ntuple = ntuple[era]
                matchs[era] = {}
                for p_nano, arr_nano in d_nano.items():
                    nano_set = set(tuple(x) for x in arr_nano)
                    for p_ntuple, arr_ntuple in d_ntuple.items():
                        ntuple_set = arr_ntuple
                        n_match = len(ntuple_set & nano_set)
                        if n_match > 0:
                            matchs[era].setdefault(p_nano, []).append(p_ntuple)
                            n_events[era][1] += n_match
                    pbar.update(1)
        print(key, str(ntuple_list_path.stem), n_events)

        
        if 1.0* n_events[key][1]/n_events[key][0]<0.8:
            temp = Path('bad.txt')
            with temp.open("a") as f: 
                f.write( ', '.join([key, str(ntuple_list_path.stem), str(n_events[key][0]), str(n_events[key][1])]) + '\n')
        else:
            temp = Path('good.txt')
            with temp.open("a") as f: 
                f.write( ', '.join([key, str(ntuple_list_path.stem), str(n_events[key][0]), str(n_events[key][1])]) + '\n')
        with open('matchs2.json', 'w') as f:
            json.dump(matchs, f)
        with open('matchs2.json', 'r') as f:
            matchs = json.load(f)

        path_matchs = {}
        for era, d in matchs.items():
            path_matchs[era] = {}
            for nano_name, ntuple_names in d.items():
                ntuple_paths = [ntuple_k_map[n] for n in ntuple_names]
                nano_path = nano_k_map[nano_name]
                path_matchs[era][nano_path] = ntuple_paths
        tmp_path = Path('merge_ntuples_commands')

        cache_path = Path('cache_2024')
        out_base = Path('/storage/af/group/phys_exotica/delayedjets/MergedNtuples/Run3/V1p19/')
        args = []
        for era in path_matchs:
            year = Path(str(era))
            dataset = ntuple_list_path.stem
            out_dir = out_base / year / dataset 
            command_out_dir = tmp_path / 'lists'/  year / dataset
            out_dir.mkdir(parents=True, exist_ok=True)
            command_out_dir.mkdir(parents=True, exist_ok=True)
            for nano_path, ntuple_paths in path_matchs[era].items():
                out_path = out_dir / nano_path.name
                inp_nano_path = command_out_dir / f'{nano_path.stem}_nano.txt'
                inp_ntuple_path = command_out_dir / f'{nano_path.stem}_ntuple.txt'
                inp_nano_path.write_text(str(nano_path).replace('root:/cmsxrootd.fnal.gov/s','root://cmsxrootd.fnal.gov//s'))
                inp_ntuple_path.write_text('\n'.join(str(x).replace('root:/cmsxrootd.fnal.gov/s','root://cmsxrootd.fnal.gov//s') for x in ntuple_paths))
                
                arg = f'{inp_ntuple_path} {inp_nano_path} {out_path} {cache_path}'
                args.append(arg)
        inp_file = tmp_path / year / Path(ntuple_list_path.stem)
        (tmp_path / year).mkdir(parents=True, exist_ok=True)
        inp_file.write_text('\n'.join(args))
        # exit()
