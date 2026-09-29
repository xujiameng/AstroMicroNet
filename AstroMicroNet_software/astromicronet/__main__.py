"""Command-line entry for matrix or processed-trace inputs (.npz/.mat)."""
from pathlib import Path
import argparse, json, sys
import numpy as np
from scipy.io import loadmat, savemat
from .pipeline import Config, fit, correlation_matrix


def tuple_key(value):
    return tuple(tuple_key(x) for x in value) if isinstance(value, list) else value


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--config', type=Path)
    args = p.parse_args()
    options = json.loads(args.config.read_text(encoding='utf-8-sig')) if args.config else {}
    cfg = Config(**options.get('pipeline', {}))
    data = (loadmat(args.input, simplify_cells=True) if args.input.suffix.lower() == '.mat'
            else dict(np.load(args.input, allow_pickle=False)))
    if 'R' in data:
        R = np.asarray(data['R'], float).copy()
        np.fill_diagonal(R, 0)
    elif 'X' in data:
        R = correlation_matrix(data['X'], options.get('definition', 'xcorr'), options.get('max_lag_frames'))
    else:
        raise ValueError('Input must contain R or processed traces X (T by N)')
    seed_key = tuple_key(options.get('seed_key', 'default'))
    gate_key = tuple_key(options['gate_key']) if 'gate_key' in options else None
    result = fit(R, cfg, data.get('roi_ids'), seed_key, gate_key)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    arrays = {k: v for k, v in result.items() if k not in ['search_logs', 'config']}
    arrays['R'] = R
    if args.output.suffix.lower() == '.mat':
        savemat(args.output, arrays, do_compression=True)
    else:
        np.savez_compressed(args.output, **arrays)
    summary = dict(config=result['config'], input=str(args.input.resolve()),
                   seed_key=seed_key, gate_key=gate_key, N=len(R), C=int(result['final'].sum()),
                   gate_p=result['gate_p'], gate_accept=result['gate_accept'],
                   theta=result['theta'], all_converged=not any(x['at_budget'] for x in result['search_logs']))
    args.output.with_suffix('.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == '__main__':
    main()
