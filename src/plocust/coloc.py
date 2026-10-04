"""Colocalization of two loci from their imprints (GWAS summary statistics) and a reference LD panel.

Steps for a pair of passports on a common build:
  1. move each imprint's SNPs to the build (local offset from the anchored lead; strand-aware),
  2. match them to reference-panel SNPs (position +- tolerance, alleles must agree),
  3. impute each study's z onto all panel SNPs of the region with panel LD
     (z_t = S_to (S_oo + lambda I)^-1 z_o; targets with low imputation quality are dropped),
  4. coloc ABF (Giambartolomei et al. 2014): posterior of
     H3 = two distinct causal variants, H4 = one shared causal variant.

Imputation needs signed z (studies that report effects). An unsigned study is
used on its directly matched SNPs only.
Like coloc, this assumes at most one causal variant per locus and study.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .io import Genotypes
from .ld import _standardize
from .passport import LocusPassport

_COMP = str.maketrans("ACGT", "TGCA")


@dataclass
class ColocConfig:
    p1: float = 1e-4
    p2: float = 1e-4
    p12: float = 1e-5
    prior_sd: float = 0.15  # effect sd in phenotype-sd units (coloc default for quantitative traits)
    position_tolerance: int = 20  # bp, after moving SNPs between builds
    ridge: float = 0.1
    min_info: float = 0.6  # minimum imputation quality kept
    min_maf: float = 0.01
    default_n: int = 1000
    h4: float = 0.8  # posterior needed for a call


def imprint_on_build(p: LocusPassport, build: str) -> Optional[pd.DataFrame]:
    """The imprint's SNPs with positions and alleles expressed on `build`."""
    imp, pl = p.imprint, p.placement(build)
    if imp is None or pl is None or pl.lead_pos is None:
        return None
    src = next((x for x in p.placements if x.method == "source"), None)
    df = pd.DataFrame({"id": imp.ids, "pos": imp.pos, "ref": imp.ref, "alt": imp.alt, "z": imp.z})
    if imp.territory_start is not None:  # only this locus's own SNPs, not a neighbouring locus's peak
        df = df[(df["pos"] >= imp.territory_start) & (df["pos"] <= imp.territory_end)]
    if len(df) == 0:
        return None
    if src is not None and src.build != build:
        d = df["pos"] - src.lead_pos
        if pl.strand == "-":
            df["pos"] = pl.lead_pos - d
            df["ref"] = df["ref"].str.upper().str.translate(_COMP)
            df["alt"] = df["alt"].str.upper().str.translate(_COMP)
        else:
            df["pos"] = pl.lead_pos + d
    df["chrom"] = pl.chrom
    return df


def to_panel(df: pd.DataFrame, geno: Genotypes, tol: int) -> pd.DataFrame:
    """Match SNPs to the panel; returns rows with panel `idx` and z for the panel's counted allele (A1).

    Exact positions are matched in one merge; only the rest are searched within +-tol bp.
    """
    empty = pd.DataFrame(columns=["idx", "pos", "z_a1"])
    if df is None or df.empty:
        return empty
    chrom = df["chrom"].iat[0]
    v = geno.region(chrom, int(df["pos"].min()) - tol, int(df["pos"].max()) + tol).sort_values("pos")
    if v.empty:
        return empty
    d = df.assign(ref=df["ref"].astype(str).str.upper(), alt=df["alt"].astype(str).str.upper())
    vv = v.assign(a1=v["a1"].astype(str).str.upper(), a2=v["a2"].astype(str).str.upper())

    def aligned(m: pd.DataFrame) -> pd.DataFrame:
        ok = ((m["ref"] == m["a1"]) & (m["alt"] == m["a2"])) | ((m["ref"] == m["a2"]) & (m["alt"] == m["a1"]))
        m = m[ok]
        sign = np.where(m["alt"] == m["a1"], 1.0, -1.0)
        return pd.DataFrame({"idx": m["idx"].astype(int).to_numpy(), "pos": m["pos_panel"].astype(int).to_numpy(),
                             "z_a1": m["z"].to_numpy() * sign, "row_": m["row_"].to_numpy()})

    d = d.reset_index(drop=True).assign(row_=lambda x: np.arange(len(x)))
    exact = aligned(d.merge(vv.rename(columns={"pos": "pos_panel"}), left_on="pos", right_on="pos_panel"))
    rest = d[~d["row_"].isin(exact["row_"])]
    near = []
    if tol > 0 and len(rest):
        vp = vv["pos"].to_numpy()
        for r in rest.itertuples():
            lo, hi = np.searchsorted(vp, r.pos - tol), np.searchsorted(vp, r.pos + tol, side="right")
            if lo == hi:
                continue
            c = vv.iloc[lo:hi]
            ok = ((c["a1"] == r.ref) & (c["a2"] == r.alt)) | ((c["a1"] == r.alt) & (c["a2"] == r.ref))
            c = c[ok]
            if len(c):
                k = c.iloc[int(np.argmin(np.abs(c["pos"].to_numpy() - r.pos)))]
                near.append((int(k["idx"]), int(k["pos"]), r.z * (1.0 if r.alt == k["a1"] else -1.0), r.row_))
    out = pd.concat([exact, pd.DataFrame(near, columns=["idx", "pos", "z_a1", "row_"])], ignore_index=True)
    return out.drop(columns="row_").drop_duplicates("idx").reset_index(drop=True)


def impute(z_obs: np.ndarray, g_obs: np.ndarray, g_tgt: np.ndarray, ridge: float) -> tuple[np.ndarray, np.ndarray]:
    """Summary-statistic imputation (ImpG / DIST style): z and quality (r2) for each target SNP."""
    so, st = _standardize(g_obs), _standardize(g_tgt)
    n = so.shape[0]
    s_oo = so.T @ so / n + ridge * np.eye(so.shape[1])
    s_to = st.T @ so / n
    w = np.linalg.solve(s_oo, s_to.T).T  # targets x observed
    z = w @ z_obs
    info = np.einsum("ij,ij->i", w, s_to)
    return z, np.clip(info, 0, 1)


def _log_abf(z: np.ndarray, n: float, maf: np.ndarray, prior_sd: float) -> np.ndarray:
    v = 1.0 / (2 * n * maf * (1 - maf))
    r = prior_sd**2 / (prior_sd**2 + v)
    return 0.5 * np.log1p(-r) + 0.5 * r * z**2


def _logsum(x: np.ndarray) -> float:
    m = x.max()
    return float(m + np.log(np.exp(x - m).sum()))


def coloc_abf(z1: np.ndarray, z2: np.ndarray, n1: float, n2: float, maf: np.ndarray, cfg: ColocConfig) -> dict:
    l1, l2 = _log_abf(z1, n1, maf, cfg.prior_sd), _log_abf(z2, n2, maf, cfg.prior_sd)
    s1, s2, s12 = _logsum(l1), _logsum(l2), _logsum(l1 + l2)
    lh = np.array([
        0.0,
        np.log(cfg.p1) + s1,
        np.log(cfg.p2) + s2,
        np.log(cfg.p1) + np.log(cfg.p2) + s1 + s2 + np.log1p(-min(np.exp(s12 - s1 - s2), 1 - 1e-12)),
        np.log(cfg.p12) + s12,
    ])
    pp = np.exp(lh - _logsum(lh))
    return {f"PP.H{i}": round(float(x), 4) for i, x in enumerate(pp)}


def coloc_pair(a: LocusPassport, b: LocusPassport, build: str, geno: Genotypes,
               cfg: Optional[ColocConfig] = None) -> Optional[dict]:
    """Colocalization posterior for two passports on `build`, or None when it cannot be computed."""
    cfg = cfg or ColocConfig()
    da, db = imprint_on_build(a, build), imprint_on_build(b, build)
    if da is None or db is None or da["chrom"].iat[0] != db["chrom"].iat[0]:
        return None
    ma, mb = to_panel(da, geno, cfg.position_tolerance), to_panel(db, geno, cfg.position_tolerance)
    if len(ma) < 2 or len(mb) < 2:
        return None
    lo, hi = max(da["pos"].min(), db["pos"].min()), min(da["pos"].max(), db["pos"].max())
    if hi <= lo:
        return None
    region = geno.region(da["chrom"].iat[0], int(lo), int(hi))
    g_all = geno.read(region["idx"].to_numpy())
    with np.errstate(invalid="ignore"):
        af = np.nanmean(g_all, axis=0) / 2
    maf = np.minimum(af, 1 - af)
    keep = np.isfinite(maf) & (maf >= cfg.min_maf)
    region, g_all, maf = region[keep], g_all[:, keep], maf[keep]
    pos_index = {int(i): k for k, i in enumerate(region["idx"].to_numpy())}

    zs = []
    for p, m in ((a, ma), (b, mb)):
        z = np.full(len(region), np.nan)
        obs = m[m["idx"].isin(pos_index)]
        cols = np.array([pos_index[int(i)] for i in obs["idx"]], dtype=int)
        z[cols] = obs["z_a1"].to_numpy()
        if p.imprint.z_signed and len(cols) >= 2:
            tgt = np.setdiff1d(np.arange(len(region)), cols)
            if len(tgt):
                zi, info = impute(obs["z_a1"].to_numpy(), g_all[:, cols], g_all[:, tgt], cfg.ridge)
                good = info >= cfg.min_info
                z[tgt[good]] = zi[good] / np.sqrt(np.maximum(info[good], 1e-6))  # rescale to unit variance
        zs.append(z)
    both = np.isfinite(zs[0]) & np.isfinite(zs[1])
    if both.sum() < 2:
        return None
    n1 = a.source.n_samples or cfg.default_n
    n2 = b.source.n_samples or cfg.default_n
    res = coloc_abf(np.abs(zs[0][both]), np.abs(zs[1][both]), n1, n2, maf[both], cfg)
    res["coloc_snps"] = int(both.sum())
    res["coloc_call"] = ("same" if res["PP.H4"] >= cfg.h4 else "distinct" if res["PP.H3"] >= cfg.h4
                         else "inconclusive")
    return res
