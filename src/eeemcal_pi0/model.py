# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/model.py and Axel Perez Ruiz's master_thesis
# pi0_calibration.ipynb, commit c36e96ce463b3bfc29abac40243cf5825400c5d1,
# cell 3 minimize_eps and cell 4 iteration. The negative sign of D was already upstream.
# Mathematical method: Jones et al., DOI 10.1016/j.nima.2006.07.061, Section 3.
# Attribution and source history: THIRD_PARTY_NOTICES.md.
"""Fixed-shower pi0 calibration; energies GeV, positions mm, cell IDs uint64.

The constrained normal-equation kernel is adapted from Axel Perez's
pi0_calibration.ipynb (c36e96ce), following Jones et al., NIM A566 (2006),
Eqs. (1)--(9). See ../../NOTICE.md for corrections and attribution.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse

M0 = 0.1349768


@dataclass
class Sample:
    arrays: dict
    showers: tuple

    @property
    def n(self):
        return len(self.arrays['energy'])

    @property
    def k(self):
        return len(self.arrays['cell_id'])

    def save(self, path):
        values = dict(self.arrays)
        for j, c in enumerate(self.showers):
            c = c.tocsr()
            values.update({f'c{j}_data': c.data, f'c{j}_indices': c.indices,
                           f'c{j}_indptr': c.indptr})
        np.savez_compressed(path, **values)

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as f:
            a = {k: f[k] for k in f.files}
        shape = (len(a['energy']), len(a['cell_id']))
        c = tuple(sparse.csr_matrix((a.pop(f'c{j}_data'), a.pop(f'c{j}_indices'),
                                    a.pop(f'c{j}_indptr')), shape=shape) for j in range(2))
        return cls(a, c)


def angular_factor(position, vertex):
    d = np.asarray(position, dtype=float) - np.asarray(vertex, dtype=float)[:, None, :]
    d /= np.linalg.norm(d, axis=2)[:, :, None]
    # Squared separation avoids cancellation for small opening angles.
    return np.sum((d[:, 0] - d[:, 1])**2, axis=1)  # 2*(1-cos(theta))


def evaluate(sample, gains, rows=None, jacobian=False):
    """All shower energy, including frozen channels, stays in denominators."""
    rows = slice(None) if rows is None else rows
    original = [c[rows] for c in sample.showers]
    current = [c.multiply(gains).tocsr() for c in original]
    sums = np.column_stack([np.asarray(c.sum(axis=1)).ravel() for c in current])
    sums0 = np.column_stack([np.asarray(c.sum(axis=1)).ravel() for c in original])
    if np.any(sums <= 0) or np.any(sums0 <= 0):
        raise ValueError('Each shower must have positive assigned energy')
    p = sample.arrays['energy'][rows] * sums / sums0
    m2 = p[:, 0] * p[:, 1] * sample.arrays['angle_factor'][rows]
    if not jacobian:
        return p, m2
    j = sum(c.multiply((m2 / sums[:, q])[:, None]) for q, c in enumerate(current))
    return p, m2, j.tocsr()


def equations(sample, gains, rows, active, chunk=2048, target=M0):
    """Merge additive C,D,L,B, never independently fitted shard gain maps."""
    k = len(active)
    C, D, L, B = np.zeros((k, k)), np.zeros(k), np.zeros(k), 0.0
    for start in range(0, len(rows), chunk):
        _, m2, J = evaluate(sample, gains, rows[start:start+chunk], jacobian=True)
        J = J[:, active]
        r = m2 - target**2
        C += (J.T @ J).toarray()
        D -= np.asarray(J.T @ r).ravel()
        L += np.asarray(J.sum(axis=0)).ravel()
        B += float(r.sum())
    return C, D, L, B


def constrained_step(C, D, L, B, rcond):
    """Unweighted paper solve; eigenvalue cutoff is not an uncertainty estimate."""
    eig, Q = np.linalg.eigh(C)
    keep = (eig > 0) & (eig > rcond * eig[-1])
    if not np.any(keep):
        raise ValueError('No supported eigenmodes')
    V, e = Q[:, keep], eig[keep]
    invD, invL = V @ ((V.T @ D) / e), V @ ((V.T @ L) / e)
    denom = float(L @ invL)
    if not np.isfinite(denom) or denom <= np.finfo(float).tiny:
        raise ValueError('Mass constraint has no retained support')
    lam = (B + L @ invD) / denom
    eps = invD - lam * invL
    return eps, {'retained_modes': int(keep.sum()), 'candidate_modes': len(eig),
                 'eigen_min_retained': float(e[0]), 'eigen_max': float(e[-1]),
                 'lambda': float(lam),
                 'linear_constraint_residual': float(B + L @ eps)}


def usage_counts(sample, rows):
    seen = (sample.showers[0][rows] + sample.showers[1][rows]).tocsr()
    seen.data = (seen.data > 0).astype(float)
    return np.asarray(seen.sum(axis=0)).ravel().astype(int)


def metrics(m2, target=M0):
    if not len(m2):
        return {'n': 0}
    mass = np.sqrt(np.maximum(m2, 0))
    return {'n': len(m2), 'mass_mean_GeV': float(mass.mean()),
            'mass_std_GeV': float(mass.std()), 'mass_median_GeV': float(np.median(mass)),
            'mass_q16_q84_GeV': np.quantile(mass, [.16, .84]).tolist(),
            'm2_bias_GeV2': float(np.mean(m2-target**2)),
            'm2_residual_std_GeV2': float(np.std(m2-target**2)),
            'm2_mse_GeV4': float(np.mean((m2-target**2)**2))}


def fit(sample, rows, min_usage=10, rcond=1e-3, iterations=8, max_step=.15,
        step_tolerance=1e-4, bias_tolerance=1e-6, chunk=2048, target=M0,
        active_indices=None):
    gains = np.ones(sample.k)
    usage = usage_counts(sample, rows)
    # Statistics studies may supply an eligibility set counted from independent
    # events and meaningful shower fractions. All columns remain in evaluate.
    active = np.flatnonzero(usage >= min_usage) if active_indices is None else np.asarray(active_indices, dtype=int)
    if len(np.unique(active)) != len(active) or np.any(active < 0) or np.any(active >= sample.k):
        raise ValueError('Invalid active channel indices')
    history = [dict(iteration=0, max_gain_step=0., **metrics(evaluate(sample, gains, rows)[1], target))]
    if len(rows) == 0 or len(active) == 0:
        return gains, usage, active, history, 'no_supported_channels'
    reason = 'iteration_limit'
    penalty = 0.
    for it in range(1, iterations+1):
        eq = equations(sample, gains, rows, active, chunk, target)
        eps, info = constrained_step(*eq, rcond)
        if not np.all(np.isfinite(eps)):
            raise ValueError('Nonfinite gain step')
        alpha = min(1., max_step/max(float(np.max(np.abs(eps))), 1e-30))
        # The equality-constrained solution need not reduce unconstrained SSE
        # when some energy is frozen. Use an exact L1 constraint merit instead:
        # phi = .5*mean(r**2) + rho*abs(mean(r)), rho > abs(lambda).
        # With L@eps=-B its directional derivative is negative away from a
        # stationary feasible point. This changes step acceptance, not C/D/L/B.
        penalty = max(penalty, 1.1*abs(info['lambda']) + 1e-12)
        old = .5*history[-1]['m2_mse_GeV4'] + penalty*abs(history[-1]['m2_bias_GeV2'])
        slope = (-float(eq[1] @ eps) - penalty*abs(eq[3]))/len(rows)
        for _ in range(25):
            trial = gains.copy()
            trial[active] *= 1 + alpha*eps
            if np.all(trial > 0) and np.all(np.isfinite(trial)):
                m = metrics(evaluate(sample, trial, rows)[1], target)
                merit = .5*m['m2_mse_GeV4'] + penalty*abs(m['m2_bias_GeV2'])
                if merit <= old + 1e-4*alpha*slope + max(1e-28, abs(old)*1e-12):
                    break
            alpha *= .5
        else:
            reason = 'line_search_stalled'
            break
        step = float(np.max(np.abs(alpha*eps)))
        spread_change = abs(m['m2_residual_std_GeV2'] - history[-1]['m2_residual_std_GeV2'])
        gains = trial
        history.append(dict(iteration=it, alpha=alpha, max_gain_step=step,
                            cumulative_gain_min=float(gains.min()),
                            cumulative_gain_max=float(gains.max()),
                            merit_GeV4=merit, merit_weight_GeV2=penalty,
                            merit_slope_GeV4=slope, **info, **m))
        if step < step_tolerance and abs(m['m2_bias_GeV2']) < bias_tolerance and spread_change < bias_tolerance:
            reason = 'converged'
            break
    return gains, usage, active, history, reason
