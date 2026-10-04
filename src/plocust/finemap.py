"""Fine-mapping from summary statistics: how many signals, and which variants.

`susie_rss` is the Sum of Single Effects model (Wang et al. 2020) in its summary-statistics form
(Zou et al. 2022): z-scores plus an LD matrix R, residual variance fixed at 1. It splits a region into up to
L single-effect components; each kept component is one signal with a credible set.

`susie_rss_multi` fits one model to several populations at once (e.g. indica and japonica panels): a signal's
causal variant is shared, its effect size may differ, and each population brings its own z and LD. Different LD
patterns break ties between variants that are linked in one population only, which narrows credible sets
(the idea behind multi-ancestry fine-mapping such as SuSiEx / MESuSiE).

Both accept prior weights per variant, e.g. from functional annotation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np


@dataclass
class Signal:
    variants: np.ndarray  # indices of the credible set, most probable first
    pip: np.ndarray  # posterior probability of each credible-set variant for this signal
    coverage: float
    purity: float  # min |r| between credible-set variants (in the first population's LD)
    log_bf: float


@dataclass
class FineMap:
    signals: list[Signal]
    pip: np.ndarray  # overall posterior inclusion probability per variant
    alpha: np.ndarray  # L x p
    prior_variance: np.ndarray  # L (x populations)
    converged: bool
    iterations: int
    notes: list[str] = field(default_factory=list)


def _norm_logpdf_ratio(b, s2, v):
    """log N(b; 0, v + s2) - log N(b; 0, s2)."""
    return 0.5 * np.log(s2 / (v + s2)) + 0.5 * b * b * (1.0 / s2 - 1.0 / (v + s2))


def _estimate_v(b, s2, log_w, grid):
    """Prior variance maximising the single-effect marginal likelihood (grid search, 0 allowed)."""
    best_v, best_ll = 0.0, 0.0  # V = 0: no effect, log-likelihood ratio 0
    for v in grid:
        lbf = _norm_logpdf_ratio(b, s2, v)
        m = (lbf + log_w).max()
        ll = m + np.log(np.exp(lbf + log_w - m).sum())
        if ll > best_ll:
            best_v, best_ll = v, ll
    return best_v


def susie_rss_multi(zs: Sequence[np.ndarray], Rs: Sequence[np.ndarray], ns: Sequence[int], L: int = 10,
                    prior_weights: Optional[np.ndarray] = None, coverage: float = 0.95, min_purity: float = 0.5,
                    max_iter: int = 100, tol: float = 1e-4) -> FineMap:
    """SuSiE on several populations with shared causal variants. zs[k], Rs[k] for population k."""
    K = len(zs)
    p = len(zs[0])
    w = np.full(p, 1.0 / p) if prior_weights is None else np.asarray(prior_weights, float) / np.sum(prior_weights)
    log_w = np.log(np.maximum(w, 1e-300))
    # sufficient statistics on the standardised scale: XtX = (n-1) R, Xty = sqrt(n-1) z, sigma2 = 1
    XtX = [np.asarray(R, float) * (n - 1) for R, n in zip(Rs, ns)]
    Xty = [np.asarray(z, float) * np.sqrt(n - 1) for z, n in zip(zs, ns)]
    d = [np.diag(A).copy() for A in XtX]
    for k in range(K):
        d[k][d[k] <= 0] = 1e-8
    alpha = np.full((L, p), 1.0 / p)
    mu = np.zeros((K, L, p))
    s2 = np.zeros((K, L, p))
    V = np.zeros((K, L))
    grid = np.concatenate([[1e-4, 1e-3, 3e-3], np.geomspace(0.01, 1.0, 12)])
    fitted = [XtX[k] @ (alpha * mu[k]).sum(0) for k in range(K)]
    converged, it = False, 0
    for it in range(1, max_iter + 1):
        old = alpha.copy()
        for l in range(L):
            lbf = np.zeros(p)
            for k in range(K):
                fitted[k] -= XtX[k] @ (alpha[l] * mu[k, l])
                r = Xty[k] - fitted[k]
                b, shat2 = r / d[k], 1.0 / d[k]
                V[k, l] = _estimate_v(b, shat2, log_w, grid)
                v = max(V[k, l], 1e-12)
                lbf += _norm_logpdf_ratio(b, shat2, v) if V[k, l] > 0 else 0.0
                s2[k, l] = 1.0 / (1.0 / v + d[k])
                mu[k, l] = s2[k, l] * r
            x = lbf + log_w
            alpha[l] = np.exp(x - x.max())
            alpha[l] /= alpha[l].sum()
            if np.all(V[:, l] == 0):
                alpha[l] = w.copy()
                mu[:, l] = 0.0
            for k in range(K):
                fitted[k] += XtX[k] @ (alpha[l] * mu[k, l])
        if np.abs(alpha - old).max() < tol:
            converged = True
            break

    signals = []
    R0 = np.asarray(Rs[0], float)
    seen = set()
    for l in range(L):
        if np.all(V[:, l] == 0):
            continue
        order = np.argsort(-alpha[l])
        k = int(np.searchsorted(np.cumsum(alpha[l][order]), coverage) + 1)
        cs = order[: min(k, p)]
        key = tuple(sorted(cs.tolist()))
        if key in seen:
            continue
        purity = float(np.abs(R0[np.ix_(cs, cs)]).min()) if len(cs) > 1 else 1.0
        if purity < min_purity:
            continue
        seen.add(key)
        lbf_l = sum(_norm_logpdf_ratio(Xty[j] / d[j] - (fitted[j] - XtX[j] @ (alpha[l] * mu[j, l])) / d[j], 1.0 / d[j],
                                       max(V[j, l], 1e-12)) for j in range(K))
        x = lbf_l + log_w
        signals.append(Signal(variants=cs, pip=alpha[l][cs], coverage=float(alpha[l][cs].sum()), purity=purity,
                              log_bf=float(x.max() + np.log(np.exp(x - x.max()).sum()))))
    pip = 1.0 - np.prod(1.0 - alpha[[l for l in range(L) if not np.all(V[:, l] == 0)]], axis=0) if signals else np.zeros(p)
    signals.sort(key=lambda s: -s.log_bf)
    notes = [] if converged else [f"did not converge in {max_iter} iterations"]
    return FineMap(signals, pip, alpha, V, converged, it, notes)


def susie_rss(z: np.ndarray, R: np.ndarray, n: int, **kw) -> FineMap:
    """Single-population SuSiE from z-scores and LD."""
    return susie_rss_multi([z], [R], [n], **kw)


def ld_consistency(z: np.ndarray, R: np.ndarray, ridge: float = 0.01) -> np.ndarray:
    """Standardised residual of each z given all others under LD R (large |value| = z/LD mismatch).

    The leave-one-out conditional expectation of z_j is computed from the precision matrix of R; values
    above ~3 point to allele flips, wrong positions, or LD from a population unlike the study's.
    """
    P = np.linalg.inv(np.asarray(R, float) + ridge * np.eye(len(z)))
    return (P @ z) / np.sqrt(np.diag(P))
