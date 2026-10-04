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


# ------------------------------------------------------------------ regions of GWAS loci


def region_z(ss, chrom: str, start: int, end: int):
    """Signed z for the alt allele of each SNP in [start, end] (sign from beta when reported)."""
    from .identify import signed_z
    from .io import normalize_chrom

    r = ss[(ss["chrom"] == normalize_chrom(chrom)) & (ss["pos"] >= start) & (ss["pos"] <= end)]
    z, signed = signed_z(r)
    return r.assign(z=z), signed


def finemap_region(ss, chrom: str, start: int, end: int, geno, n: int, L: int = 10, max_snps: int = 2500,
                   min_maf: float = 0.01, panel_name: str = "panel", covariates: Optional[np.ndarray] = None,
                   min_abs_z: Optional[float] = 4.42, **kw):
    """Fine-map one region of one study: summary statistics `ss`, LD from genotypes `geno` (same build).

    SNPs are matched to the panel by position and alleles (z flipped to the panel's counted allele).
    Above max_snps, the most significant SNPs are kept.

    `covariates` (samples x k, same sample order as `geno`): the GWAS covariates, e.g. its PCs. Mixed-model
    and PC-adjusted z-scores describe genotypes with population structure removed, so LD is computed on
    genotypes residualised on the same covariates; raw LD makes SuSiE invent signals.
    `min_abs_z`: a signal is kept only if one of its credible-set SNPs reaches this |z| (4.42 ~ p 1e-5).
    Returns (FineMap, matched SNP table, LD matrix) or None.
    """
    from .coloc import to_panel

    r, signed = region_z(ss, chrom, start, end)
    if not signed or len(r) < 2:
        return None
    df = r.rename(columns={"z": "z"})[["chrom", "pos", "ref", "alt", "z"]]
    m = to_panel(df.assign(chrom=chrom), geno, tol=0)
    if len(m) < 2:
        return None
    g = geno.read(m["idx"].to_numpy())
    with np.errstate(invalid="ignore"):
        af = np.nanmean(g, axis=0) / 2
    keep = np.isfinite(af) & (np.minimum(af, 1 - af) >= min_maf)
    m, g = m[keep].reset_index(drop=True), g[:, keep]
    if len(m) > max_snps:
        top = np.sort(np.argsort(-np.abs(m["z_a1"].to_numpy()))[:max_snps])
        m, g = m.iloc[top].reset_index(drop=True), g[:, top]
    if len(m) < 2:
        return None
    from .ld import _standardize

    x = _standardize(g)
    if covariates is not None:
        c = np.column_stack([np.ones(len(x)), np.asarray(covariates, float)])
        x = x - c @ np.linalg.lstsq(c, x, rcond=None)[0]
        sd = x.std(0)
        x = x / np.where(sd == 0, 1, sd)
    R = (x.T @ x) / x.shape[0]
    np.fill_diagonal(R, 1.0)
    fm = susie_rss(m["z_a1"].to_numpy(), R, n=n, L=L, **kw)
    if min_abs_z is not None:
        z = np.abs(m["z_a1"].to_numpy())
        fm.signals = [s for s in fm.signals if z[s.variants].max() >= min_abs_z]
    resid = ld_consistency(m["z_a1"].to_numpy(), R)
    fm.notes.append(f"ld_outliers={int((np.abs(resid) > 4).sum())}")
    return fm, m.assign(resid=resid), R


def to_fine_mapping(fm, m, start: int, end: int, L: int, ld_source: list[str], method: str = "susie_rss",
                    lead_pos: Optional[int] = None, R: Optional[np.ndarray] = None):
    """Convert a FineMap on SNP table `m` (pos, idx) to the passport schema."""
    from .passport import CredibleVariant, FineMapping, FineMapSignal

    pos = m["pos"].to_numpy()
    signals = [FineMapSignal(log_bf=round(s.log_bf, 3), purity=round(min(max(s.purity, 0.0), 1.0), 3),
                             coverage=round(min(s.coverage, 1.0), 4),
                             variants=[CredibleVariant(pos=int(pos[v]), pip=round(float(min(pp, 1.0)), 5))
                                       for v, pp in zip(s.variants, s.pip)]) for s in fm.signals]
    lead_signal = None
    if lead_pos is not None and fm.signals:
        hit = [i for i, s in enumerate(fm.signals) if lead_pos in set(pos[s.variants].tolist())]
        if hit:
            lead_signal = hit[0]
        elif R is not None and lead_pos in set(pos.tolist()):
            j = int(np.flatnonzero(pos == lead_pos)[0])
            r2 = [float(np.max(R[j, s.variants] ** 2)) for s in fm.signals]
            lead_signal = int(np.argmax(r2)) if max(r2) >= 0.5 else None
    outliers = next((int(n.split("=")[1]) for n in fm.notes if n.startswith("ld_outliers=")), 0)
    return FineMapping(method=method, ld_source=ld_source, region_start=int(start), region_end=int(end),
                       n_snps=len(m), max_signals=L, converged=fm.converged, ld_outliers=outliers, signals=signals,
                       lead_signal=lead_signal)


def loci_regions(passports, build: str, flank: int = 250_000, max_len: int = 2_000_000):
    """Group a study's loci into regions: lead +- flank, overlapping windows merged (capped at max_len)."""
    spans = []
    for p in passports:
        pl = p.placement(build)
        if pl is not None and pl.lead_pos is not None:
            spans.append((pl.chrom, max(1, pl.lead_pos - flank), pl.lead_pos + flank, p))
    spans.sort(key=lambda t: (t[0], t[1]))
    regions = []
    for chrom, a, b, p in spans:
        if regions and regions[-1][0] == chrom and a <= regions[-1][2] and b - regions[-1][1] <= max_len:
            regions[-1][2] = max(regions[-1][2], b)
            regions[-1][3].append(p)
        else:
            regions.append([chrom, a, b, [p]])
    return [tuple(r) for r in regions]


def finemap_passports(passports, ss, geno, build: str, n: Optional[int] = None, L: int = 10,
                      covariates: Optional[np.ndarray] = None, ld_source: str = "panel", flank: int = 250_000):
    """Fine-map the regions of a study's loci and attach the result to every passport. Returns region rows."""
    rows = []
    for chrom, a, b, members in loci_regions(passports, build, flank=flank):
        nn = n or next((p.source.n_samples for p in members if p.source.n_samples), None) or 1000
        res = finemap_region(ss, chrom, a, b, geno, n=nn, L=L, covariates=covariates)
        if res is None:
            continue
        fm, m, R = res
        for p in members:
            p.fine_mapping = to_fine_mapping(fm, m, a, b, L, [ld_source], lead_pos=p.placement(build).lead_pos, R=R)
        rows.append({"chrom": chrom, "start": a, "end": b, "loci": len(members), "snps": len(m),
                     "signals": len(fm.signals), "cs_sizes": ",".join(str(len(s.variants)) for s in fm.signals),
                     "ld_outliers": int((np.abs(m["resid"]) > 4).sum())})
    return rows
