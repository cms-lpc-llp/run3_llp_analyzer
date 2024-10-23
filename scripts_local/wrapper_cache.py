from __future__ import annotations
from multiprocessing import Pool
from pathlib import Path

import click
import numpy as np
import uproot as up
import h5py as h5
import os
base_path = '/storage/af/user/christiw/login-1/christiw/LLP/Run3/CMSSW_10_6_30'

##### run ntupler ######
list_path = Path(base_path + '/src/run3_llp_analyzer/lists/displacedJetMuonNtuple/V1p19/')
campaigns = {
    # 'Data2022': 'EXOCSCCluster'
    #'Data2023': 'EXOCSCCluster',
    'Data2024': 'Muon',
    # 'MC_Summer22': 'ggH',
    # 'MC_Summer22EE': 'ggH',
    # 'MC_Summer23'  : 'ggH',
    # 'MC_Summer23BPix': 'ggH',
} 
for k, v in campaigns.items():
    paths = list(list_path.glob(f"{k}/**/*{v}*"))
    for p in paths:
        print(f"running cache for: {p}")
        cache_path = str(p.parents[0]).replace('lists','/data/cache/')
        cache_name = p.stem
        print(cache_path, cache_name)
        Path(cache_path).mkdir(exist_ok=True, parents=True)
        os.system(f"python3 cache.py -i {p} -o {cache_path}/{cache_name} -j 4")
        os.system(f"python3 extract_cache.py -i {cache_path}/{cache_name} -o cache_root")

##### run nanoAOD ######
list_path = Path(base_path + '/src/run3_llp_analyzer/lists/nanoAOD/')
campaigns = {
    # '2022': 'DisplacedJet',
    # '2023': 'Muon',
    #'2024':'Muon',
    #'MC_Summer22': 'ggH',
    #'MC_Summer22EE': 'ggH',
    #'MC_Summer23'  : 'ggH',
    #'MC_Summer23BPix': 'ggH',
} 

for k, v in campaigns.items():
    paths = list(list_path.glob(f"{k}/**/*{v}*"))
    for p in paths:
        print(f"running cache for: {p}")
        cache_path = str(p.parents[0]).replace('lists','/data/cache/')
        cache_name = p.stem
        Path(cache_path).mkdir(exist_ok=True, parents=True)
        print(f"python3 cache.py -i {p} -o {cache_path}/{cache_name} -j 12")
        os.system(f"python3 cache.py -i {p} -o {cache_path}/{cache_name} -j 12")
        os.system(f"python3 extract_cache.py -i {cache_path}/{cache_name} -o cache_root")
