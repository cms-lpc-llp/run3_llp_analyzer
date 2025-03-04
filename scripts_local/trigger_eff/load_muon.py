import argparse
from pathlib import Path

import awkward as ak
import h5py as h5
import numpy as np
import uproot as up
import os
import re
import gzip
import json

m = re.compile(r'Run20\d{2}[A-Z]')

keep = [
    'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight',
    'HLT_IsoMu27',
    'PuppiMET_pt',
    'PuppiMET_phi',
    'Muon_pfRelIso04_all',
    'Muon_tightId',
    'Muon_phi',

    'Flag_goodVertices',
    'Flag_globalSuperTightHalo2016Filter',
    'Flag_EcalDeadCellTriggerPrimitiveFilter',
    'Flag_BadPFMuonFilter',
    'Flag_BadPFMuonDzFilter',
    'Flag_hfNoisyHitsFilter',
    'Flag_eeBadScFilter',
]

mc_map = {
    'Run2022B': 'Summer22',
    'Run2022C': 'Summer22',
    'Run2022D': 'Summer22',
    'Run2022E': 'Summer22EE',
    'Run2022F': 'Summer22EE',
    'Run2022G': 'Summer22EE',
    'Run2023B': 'Summer23',
    'Run2023C': 'Summer23',
    'Run2023D': 'Summer23BPix',
    'Run2024': 'Summer24',
}


def load_vet_map(veto_json_path: str):
    with gzip.open(veto_json_path, 'rt') as f:
        jet_veto_map = json.load(f)
    hist = jet_veto_map['corrections'][0]['data']['content'][0]['value']
    eta_edges, phi_edges = hist['edges']
    val = np.array(hist['content']) == 0
    eta_edges, phi_edges = np.array(eta_edges), np.array(phi_edges)
    return eta_edges, phi_edges, val

def delta_phi(phi1, phi2):
    dphi1 = np.abs(phi1 - phi2)
    dphi2 = 2 * np.pi - dphi1
    return np.minimum(dphi1, dphi2)
    
def compute_jet_mask(f, campaign: str):
    veto_json_path = f'../../data/JetVetoMap/{campaign}.json.gz'
    eta_edges, phi_edges, val = load_vet_map(veto_json_path)
    eta = f['/Events/Jet_eta'].array()
    pt = f['/Events/Jet_pt'].array()
    phi = f['/Events/Jet_phi'].array()
    met_phi = f['/Events/MET_phi'].array()
    layout = ak.num(eta)

    np_eta = np.asarray(ak.ravel(eta))
    np_phi = np.asarray(ak.ravel(phi))

    eta_bin, phi_bin = np.searchsorted(eta_edges, np_eta), np.searchsorted(phi_edges, np_phi)
    eta_bin, phi_bin = np.clip(eta_bin, 1, len(eta_edges) - 1) - 1, np.clip(phi_bin, 1, len(phi_edges) - 1) - 1
    idx = eta_bin * (len(phi_edges) - 1) + phi_bin

    jet_veto = np.asarray(ak.all(ak.unflatten(val[idx], layout), axis=-1))

    if campaign == 'Summer24':
        bad_eval_cali_mask = f['/Events/Flag_ecalBadCalibFilter'].array().to_numpy()
    else:
        run_num = f['/Events/run'].array()
        run_num = ak.broadcast_arrays(run_num, eta)[0]
        met = f['/Events/MET_pt'].array()
        met = ak.broadcast_arrays(met, eta)[0]
        neEmEF, chEmEF = f['/Events/Jet_neEmEF'].array(), f['/Events/Jet_chEmEF'].array()
        _bad_ecal_cali_reject = (
            (362433 <= run_num) & (run_num <= 367144) &
            (met > 100) &
            (pt > 50) &
            (-0.5 > eta) & (eta > -0.1) &
            (-2.1 > phi) & (phi > -1.8) &
            ((neEmEF > 0.9) | (chEmEF > 0.9)) &
            (delta_phi(met_phi, phi) > 2.9)
        )
        bad_eval_cali_reject = ak.any(_bad_ecal_cali_reject, axis=-1).to_numpy()
        bad_eval_cali_mask = ~bad_eval_cali_reject

    bad_eval_cali_mask = np.asarray(bad_eval_cali_mask)
    
    return jet_veto & bad_eval_cali_mask


def load_file(fpath: str, campaign: str, is_data: bool) -> dict[str, np.ndarray]:

    if is_data:
        keep.extend(('event', 'run', 'luminosityBlock'))
    else:
        keep.extend(('Generator_weight', 'Pileup_nPU'))

    with up.open(fpath) as f:  # type: ignore
        mu_pt = f['/Events/Muon_pt'].array()  # type: ignore
        muon_mask = (ak.num(mu_pt) == 1).to_numpy()
        jet_veto_mask = compute_jet_mask(up.open(args.input), campaign)
        mask = muon_mask & jet_veto_mask
        mu_pt: np.ndarray = mu_pt[mask].to_numpy()
        arrs = {k: f[f'/Events/{k}'].array()[mask].to_numpy().ravel() for k in keep}  # type: ignore
    arrs['Muon_pt'] = mu_pt.ravel()
    return arrs


def save(arrs: dict[str, np.ndarray], fpath: str | Path):
    with h5.File(fpath, 'w') as f:
        for k, v in arrs.items():
            f.create_dataset(k, data=v, compression='gzip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str)
    parser.add_argument('output', type=str)
    args = parser.parse_args()

    if Path(args.output).exists():
        exit(0)

    # Use local cache if exists
    fpath = args.input
    alt = Path('/tmp') / Path(fpath).name
    if alt.exists():
        fpath = str(alt)

    # Check if data or MC, get campaign
    run = m.findall(args.input)
    is_data = len(run) > 0
    if is_data:
        campaign = mc_map[run[0]]
        # print(f'Processing {run[0]} ({campaign})')
    else:
        campaign = [k for k in mc_map.values() if k in args.input][0]
        # print(f'Processing MC ({campaign}).')

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    arrs = load_file(fpath, campaign, is_data)
    save(arrs, args.output + '.tmp')
    os.rename(args.output + '.tmp', args.output)
