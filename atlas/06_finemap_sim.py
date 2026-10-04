"""Atlas step 6: does cross-subpopulation fine-mapping work on real rice genotypes?

Real 3K genotypes (dense filtered SNPs), split by PC1 into panel A (n~1,400) and panel B (n~850).
In random 150 kb regions, 1-3 causal SNPs (common in both panels) get effects; z-scores are computed in each panel.

Compared:
  A alone, B alone, A+B jointly (shared causal variants, own LD)            -> does joint fine-mapping help?
  A's z with B's LD ("published results only, LD borrowed from another panel") -> how much does LD mismatch hurt?
Metrics per method: credible sets that contain a causal SNP (coverage), credible-set size,
causal SNPs captured (power), and false signals.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from plocust.finemap import susie_rss, susie_rss_multi
from plocust.io import Genotypes

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RES = ROOT / "atlas" / "results"
N_REGIONS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
WINDOW = 150_000
H2_EACH = 0.02
RNG = np.random.default_rng(7)


def standardise(g):
    g = np.where(np.isnan(g), np.nanmean(g, axis=0), g)
    sd = g.std(0)
    return (g - g.mean(0)) / np.where(sd == 0, 1, sd), sd > 0


def z_scores(x, y):
    y = (y - y.mean()) / y.std()
    return x.T @ y / np.sqrt(len(y))


def evaluate(fm, causal, label, region):
    sets = [set(map(int, s.variants)) for s in fm.signals]
    hit = [bool(s & causal) for s in sets]
    return {"region": region, "method": label, "true_causal": len(causal), "signals": len(sets),
            "true_signals": sum(hit), "false_signals": len(sets) - sum(hit),
            "causal_captured": len(causal & set().union(*sets)) if sets else 0,
            "mean_cs_size": float(np.mean([len(s) for s in sets])) if sets else np.nan,
            "median_cs_size_true": float(np.median([len(s) for s, h in zip(sets, hit) if h])) if any(hit) else np.nan}


def main():
    dense = Genotypes.open(DATA / "3k/filt/base_filtered_v0.7")
    keep_a = pd.read_csv(DATA / "sim/keep_A.txt", sep="\t", header=None)[1]
    keep_b = pd.read_csv(DATA / "sim/keep_B.txt", sep="\t", header=None)[1]
    ga, gb = dense.subset(keep_a), dense.subset(keep_b)
    v = dense.variants
    rows, done, tries = [], 0, 0
    while done < N_REGIONS and tries < N_REGIONS * 5:
        tries += 1
        s = v.sample(1, random_state=int(RNG.integers(1e9))).iloc[0]
        reg = dense.region(s["chrom"], int(s["pos"]), int(s["pos"]) + WINDOW)
        if len(reg) < 200:
            continue
        xa, oka = standardise(ga.read(reg["idx"].to_numpy()))
        xb, okb = standardise(gb.read(reg["idx"].to_numpy()))
        fa = np.nanmean(ga.read(reg["idx"].to_numpy()), 0) / 2
        fb = np.nanmean(gb.read(reg["idx"].to_numpy()), 0) / 2
        common = oka & okb & (np.minimum(fa, 1 - fa) >= 0.05) & (np.minimum(fb, 1 - fb) >= 0.05)
        keep = np.flatnonzero(common)
        if len(keep) > 1500:  # cap matrix size
            keep = np.sort(RNG.choice(keep, 1500, replace=False))
        if len(keep) < 100:
            continue
        xa, xb = xa[:, keep], xb[:, keep]
        Ra, Rb = np.corrcoef(xa, rowvar=False), np.corrcoef(xb, rowvar=False)
        k = int(RNG.integers(1, 4))
        causal = set(RNG.choice(len(keep), k, replace=False).tolist())
        c = sorted(causal)
        effects = RNG.choice([-1, 1], k) * np.sqrt(H2_EACH)
        ya = xa[:, c] @ effects + RNG.normal(0, np.sqrt(1 - H2_EACH * k), len(xa))
        yb = xb[:, c] @ (effects * RNG.uniform(0.6, 1.4, k)) + RNG.normal(0, np.sqrt(1 - H2_EACH * k), len(xb))
        za, zb = z_scores(xa, ya), z_scores(xb, yb)
        if max(np.abs(za).max(), np.abs(zb).max()) < 5.2:  # nothing to fine-map (about p > 2e-7)
            continue
        region = f"{s['chrom']}:{int(s['pos'])}"
        fits = {
            "A alone": susie_rss(za, Ra, len(xa), L=5),
            "B alone": susie_rss(zb, Rb, len(xb), L=5),
            "A+B joint": susie_rss_multi([za, zb], [Ra, Rb], [len(xa), len(xb)], L=5),
            "A with B's LD (borrowed)": susie_rss(za, Rb, len(xa), L=5),
        }
        for label, fm in fits.items():
            rows.append(evaluate(fm, causal, label, region))
        done += 1
        print(f"region {done}/{N_REGIONS} {region}: {len(keep)} SNPs, {k} causal", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "finemap_sim_regions.tsv", sep="\t", index=False)
    s = df.groupby("method", sort=False).agg(
        regions=("region", "size"), causal=("true_causal", "sum"), captured=("causal_captured", "sum"),
        signals=("signals", "sum"), true_signals=("true_signals", "sum"), false_signals=("false_signals", "sum"),
        cs_size_median=("median_cs_size_true", "median"))
    s["power"] = (s["captured"] / s["causal"]).round(3)
    s["coverage"] = (s["true_signals"] / s["signals"]).round(3)
    s.to_csv(RES / "finemap_sim_summary.tsv", sep="\t")
    print(s.to_string())


if __name__ == "__main__":
    main()
