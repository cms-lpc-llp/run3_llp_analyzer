import h5py as h5
import numpy as np
from pathlib import Path
import uproot as up
import argparse
from glob import glob


def _find_match_idx(vec1: np.ndarray, vec2: np.ndarray):
    idx1 = np.argsort(vec1)
    svec1 = vec1[idx1]
    p_locs2 = np.searchsorted(svec1, vec2).clip(0, len(vec1) - 1)
    locs2 = np.where(svec1[p_locs2] == vec2)[0]
    locs1 = idx1[p_locs2[locs2]]
    return locs1, locs2


def encode(rle1: np.ndarray, rle2: np.ndarray):
    r1, _, e1 = rle1.T
    r2, _, e2 = rle2.T
    C = max(e1.max(), e2.max()) + 1

    vec1 = r1 * C + e1
    vec2 = r2 * C + e2
    return vec1, vec2


def find_match_idx(rle1: np.ndarray, rle2: np.ndarray):
    vec1, vec2 = encode(rle1, rle2)
    return _find_match_idx(vec1, vec2)


def gen_match_object(f1s: list[str | Path], f2s: list[str | Path], r=True):
    _rle1 = []
    _l1s = []
    _src_paths_1 = []
    for f1 in f1s:
        with h5.File(f1, 'r') as f:
            _arr = np.array(f['cache'])
            _rle1.append(_arr)
            _l1s.append(len(_arr))
            _src_paths_1.append(f.attrs['key'])
    rle1 = np.concatenate(_rle1)
    if len(rle1) == 0:
        return None

    _rle2 = []
    _l2s = []
    _src_paths_2 = []
    for f2 in f2s:
        with h5.File(f2, 'r') as f:
            _arr = np.array(f['cache'])
            _rle2.append(_arr)
            _l2s.append(len(_arr))
            _src_paths_2.append(f.attrs['key'])
    rle2 = np.concatenate(_rle2)
    if len(rle2) == 0:
        return None

    _l1s, _l2s = np.array(_l1s), np.array(_l2s)
    cl1s, cl2s = np.cumsum(_l1s), np.cumsum(_l2s)

    locs1, locs2 = find_match_idx(rle1, rle2)
    if locs1.size == 0:
        return None
    assert np.all(rle1[locs1] == rle2[locs2])
    m0 = locs1.size
    # print(f'Found {m0} matches')

    bias1 = np.searchsorted(cl1s, locs1, side='right')
    bias2 = np.searchsorted(cl2s, locs2, side='right')
    mask1 = np.unique(bias1)
    mask2 = np.unique(bias2)

    rle1 = np.concatenate([_rle1[i] for i in mask1])
    rle2 = np.concatenate([_rle2[i] for i in mask2])
    l1s = _l1s[mask1]
    l2s = _l2s[mask2]
    src_paths_1 = [_src_paths_1[i] for i in mask1]
    src_paths_2 = [_src_paths_2[i] for i in mask2]

    locs1, locs2 = find_match_idx(rle1, rle2)
    assert np.all(rle1[locs1] == rle2[locs2])
    assert m0 == locs1.size
    assert np.sum(l1s) == len(rle1) and np.sum(l2s) == len(rle2)

    return {
        'inp1': {
            'paths': src_paths_1,
            'lengths': l1s,
        },
        'inp2': {
            'paths': src_paths_2,
            'lengths': l2s,
        },
        'idx': {
            '1': locs1,
            '2': locs2,
        },
    }


def save_match_object(match_obj: dict, f: str | Path):
    Path(f).parent.mkdir(parents=True, exist_ok=True)
    with up.recreate(f) as f:  # type: ignore
        for k, v in match_obj.items():
            f[k] = v  # type: ignore


def main(fs1: list[str | Path], fs2: list[str | Path], out: str | Path):
    match_obj = gen_match_object(fs1, fs2)
    if match_obj is None:
        # print('No match found')
        return 0
    save_match_object(match_obj, out)
    return len(match_obj['idx']['1'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', type=str, nargs='+', required=False, default=None)
    parser.add_argument('-i2', '--input2', type=str, nargs='+', required=False, default=None)
    parser.add_argument('-i1r', '--input1r', type=str, required=False)
    parser.add_argument('-i2r', '--input2r', type=str, required=False)
    parser.add_argument('-o', '--output', type=str, required=True)
    args = parser.parse_args()

    input1, input2 = args.input1, args.input2
    if not input1:
        input1 = glob(args.input1r + '/*')
    if not input2:
        input2 = glob(args.input2r + '/*')

    main(input1, input2, args.output)
