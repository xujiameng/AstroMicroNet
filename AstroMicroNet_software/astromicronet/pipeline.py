"""One numerical entry used by revision analyses and the MATLAB bridge.

The accepted 20260914 search/refinement kernels are copied verbatim. This
module consolidates validation, fusion and the mandatory homogeneity screen.
The screen tests exchangeability of weights, not every possible non-CP model.
"""
from dataclasses import dataclass, asdict
import hashlib
import numpy as np
from scipy.fft import rfft, irfft, next_fast_len
from scipy.stats import rankdata
from .robust_solver import solve, refine, threshold2, objective
from .stable_fusion import stable_score_order


@dataclass(frozen=True)
class Config:
    rho: float = .20
    tau: float = .20
    random_starts: int = 1024
    row_starts: int = 16
    theta_offset: float = 0.


def validate_matrix(R):
    A = np.asarray(R, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1] or len(A) < 3:
        raise ValueError('R must be a square matrix with at least three ROIs')
    if not np.isfinite(A).all() or not np.allclose(A, A.T, atol=1e-12, rtol=0):
        raise ValueError('R must be finite and symmetric')
    if not np.allclose(np.diag(A), 0, atol=1e-12, rtol=0):
        raise ValueError('R diagonal must be zero; exclude self-connections')
    return A


def screen(R, key):
    """Accepted 499-permutation screen, with private reproducible RNG."""
    seed = int.from_bytes(hashlib.sha256(str(key).encode()).digest()[:4], 'little')
    rng = np.random.default_rng(seed)
    I, J = np.triu_indices(len(R), 1)
    weights = R[I, J]
    observed = float(np.var(R.sum(1)))
    null = np.empty(499)
    for b in range(499):
        v = rng.permutation(weights)
        strengths = np.bincount(I, weights=v, minlength=len(R)) + np.bincount(J, weights=v, minlength=len(R))
        null[b] = np.var(strengths)
    p = (1 + int(np.count_nonzero(null >= observed - 1e-12))) / 500
    return dict(gate_p=p, gate_accept=bool(p <= .05), gate_observed=observed,
                gate_null_mean=float(null.mean()), gate_null_95=float(np.quantile(null, .95)))


class ScreenCache:
    """Cache only exact matrix/key pairs; callers cannot supply gate outcomes."""
    def __init__(self):
        self._values = {}

    def evaluate(self, A, key):
        digest = hashlib.sha256(np.ascontiguousarray(A).tobytes()).hexdigest()
        identity = (A.shape, digest, repr(key))
        if identity not in self._values:
            self._values[identity] = screen(A, key)
        return dict(self._values[identity])


def fit(R, config=None, roi_ids=None, seed_key='default', gate_key=None,
        candidate_pool=None, screen_cache=None):
    """Return final labels. No public option bypasses the null screen.

    candidate_pool is a reusable set of converged candidates, not final labels.
    Its shape, fixed-k constraint and original objective values are checked.
    """
    cfg = config or Config()
    if not 0 < cfg.rho < 1 or not 0 < cfg.tau <= 1:
        raise ValueError('rho must be in (0,1); tau in (0,1]')
    if cfg.random_starts < 1 or cfg.row_starts < 0:
        raise ValueError('Invalid number of starts')
    A = validate_matrix(R)
    N = len(A)
    k = int(np.floor(cfg.rho * N))
    if not 0 < k < N:
        raise ValueError('rho and N must produce 0 < floor(rho*N) < N')
    ids = np.arange(N) if roi_ids is None else np.asarray(roi_ids).reshape(-1)
    if len(ids) != N or not np.all(ids[1:] > ids[:-1]):
        raise ValueError('roi_ids must be strictly increasing immutable ROI IDs')
    theta = threshold2(A) + cfg.theta_offset
    key = ('astromicronet_gate', seed_key) if gate_key is None else gate_key
    gate = screen(A, key) if screen_cache is None else screen_cache.evaluate(A, key)
    logs = []
    if candidate_pool is None:
        solved = solve(A, cfg.rho, seed_key, cfg.random_starts, cfg.row_starts)
        C, F, logs = solved['all_candidates'], solved['all_fitness'], solved['logs']
        if any(row['at_budget'] for row in logs):
            raise RuntimeError('Search reached its swap budget before convergence')
    else:
        C, F = np.asarray(candidate_pool[0], bool), np.asarray(candidate_pool[1], float)
    expected_starts = 1 + min(cfg.row_starts, N) + cfg.random_starts
    if C.shape != (expected_starts, N) or F.shape != (expected_starts,) or not np.all(C.sum(1) == k):
        raise ValueError('Candidate pool does not match the declared fixed-k search')
    unique, inverse = np.unique(C, axis=0, return_inverse=True)
    original_fitness = np.array([objective(A, c) for c in unique])[inverse]
    if not np.isfinite(F).all() or not np.allclose(F, original_fitness, atol=1e-9, rtol=0):
        raise ValueError('Candidate fitness does not match the input matrix')
    refined = np.asarray([refine(A, c, theta, tie_policy='first')[0] for c in unique])
    assert np.all(refined <= unique)
    weights = np.maximum(F, 0)
    support = (np.average(refined[inverse], axis=0, weights=weights)
               if weights.sum() > 0 else np.zeros(N))
    order = stable_score_order(support, ids)
    eligible = support >= cfg.tau
    chosen = order[eligible[order]][:k]
    before = np.zeros(N, bool)
    before[chosen] = True
    final = before.copy() if gate['gate_accept'] else np.zeros(N, bool)
    return dict(final=final, before_gate=before, support=support, theta=theta, k=k,
                all_candidates=C, all_fitness=F, unique_candidates=unique,
                unique_refined=refined, candidate_to_unique=inverse, roi_ids=ids,
                eligibility_count_before_cap=int(eligible.sum()),
                cap_applied=bool(eligible.sum() > k),
                removed_by_cap=int(eligible.sum() - before.sum()),
                search_logs=logs, config=asdict(cfg), **gate)


def strict_unique_records(raw, centers, boundaries, df=None):
    """Only collapse exact raw-trace AND exact-geometry duplicate records."""
    raw = np.asarray(raw)
    parent = np.arange(raw.shape[1])
    for i in range(raw.shape[1]):
        for j in range(i + 1, raw.shape[1]):
            if (np.array_equal(raw[:, i], raw[:, j]) and
                np.array_equal(centers[i], centers[j]) and
                np.array_equal(boundaries[i], boundaries[j])):
                if df is not None and not np.array_equal(df[:, i], df[:, j]):
                    raise ValueError('Same raw ROI record has inconsistent processed traces')
                parent[j] = parent[i]
    keep = np.unique(parent)
    return keep, np.searchsorted(keep, parent)


def correlation_matrix(X, definition='xcorr', max_lag_frames=None):
    """T-by-N traces. xcorr matches MATLAB coeff: no demean, no abs.

    Pearson and Spearman are zero-lag signed coefficients. They are explicitly
    different definitions, not aliases of non-demeaned xcorr.
    """
    X = np.asarray(X, float)
    if X.ndim != 2 or min(X.shape) < 2 or not np.isfinite(X).all():
        raise ValueError('Finite T-by-N traces are required')
    T, N = X.shape
    if definition in ['pearson', 'spearman']:
        Y = rankdata(X, axis=0) if definition == 'spearman' else X.copy()
        Y = Y - Y.mean(axis=0)
        norms = np.linalg.norm(Y, axis=0)
        if np.any(norms == 0):
            raise ValueError('Constant traces cannot define Pearson/Spearman edges')
        Y /= norms
        R = Y.T @ Y
    elif definition == 'xcorr':
        norms = np.linalg.norm(X, axis=0)
        if np.any(norms == 0):
            raise ValueError('Zero-energy traces cannot define xcorr coeff')
        lag = T - 1 if max_lag_frames is None else int(max_lag_frames)
        if not 0 <= lag < T:
            raise ValueError('Lag must lie in [0,T-1]')
        Y = X / norms
        if lag == 0:
            R = Y.T @ Y
        else:
            nfft = next_fast_len(2 * T - 1)
            Z = rfft(Y, n=nfft, axis=0)
            R = np.zeros((N, N))
            indices = np.r_[np.arange(lag + 1), np.arange(nfft - lag, nfft)]
            for i in range(N - 1):
                corr = irfft(Z[:, i, None] * Z[:, i+1:].conj(), n=nfft, axis=0)
                values = corr[indices].max(axis=0)
                R[i, i+1:] = values
                R[i+1:, i] = values
    else:
        raise ValueError('Unknown connectivity definition: ' + definition)
    np.fill_diagonal(R, 0)
    return validate_matrix(R)


def categories(A, core, theta):
    I, J = np.triu_indices(len(A), 1)
    kind = core[I].astype(int) + core[J].astype(int)
    out = dict(N=len(A), C=int(core.sum()), core_fraction=float(core.mean()))
    for code, name in [(2, 'CC'), (1, 'CP'), (0, 'PP')]:
        v = A[I[kind == code], J[kind == code]]
        retained = v[v >= theta]
        out[name + '_density'] = float(np.mean(v >= theta)) if len(v) else np.nan
        out[name + '_fc_allpairs'] = float(v.mean()) if len(v) else np.nan
        out[name + '_fc_retained'] = float(retained.mean()) if len(retained) else np.nan
    return out


def overlap(reference, current):
    a, b = np.asarray(reference, bool), np.asarray(current, bool)
    inter, union = int((a & b).sum()), int((a | b).sum())
    return dict(jaccard=inter/union if union else 1., intersection=inter,
                n_default=int(a.sum()), n_variant=int(b.sum()),
                default_retention=inter/a.sum() if a.any() else np.nan,
                variant_retention=inter/b.sum() if b.any() else np.nan)
