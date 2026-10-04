"""Atlas step 7: fine-map every GWAS locus (SuSiE from summary statistics + in-sample LD).

Loci of one trait are grouped into regions (lead +- 250 kb, overlaps merged, <= 2 Mb). Each region is fine-mapped
once with the GWAS panel's own genotypes as LD (3K: core SNPs on IRGSP-1.0; RDP1: 44K on MSU6). Every passport
gets the fine-mapping of its region and the index of the signal its lead belongs to.
Output: passports in data/atlas/passports_fm/, summary tables in atlas/results/.
"""

import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from plocust.finemap import finemap_region, loci_regions, to_fine_mapping
from plocust.io import Genotypes, read_sumstats
from plocust.passport import read_passports, write_passports

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PASS = DATA / "atlas" / "passports"
OUT = DATA / "atlas" / "passports_fm"
RES = ROOT / "atlas" / "results"
L = 10

PANELS = {"3K": dict(build="IRGSP-1.0", geno="3k/core3k", ld="3K core SNPs (365K), PC-adjusted", pcs="3k/pca3k.eigenvec", k=5),
          "RDP1": dict(build="MSU6", geno="rdp1/rdp1", ld="RDP1 44K array, PC-adjusted", pcs="rdp1/rpca.eigenvec", k=4)}


def panel_with_pcs(cfg):
    """The GWAS samples of the panel, and their PCs in the same order (the GWAS covariates)."""
    pcs = pd.read_csv(DATA / cfg["pcs"], sep="\t").rename(columns={"#FID": "FID"})
    geno = Genotypes.open(DATA / cfg["geno"]).subset(pcs["IID"].astype(str))
    order = geno.samples[geno._keep]["iid"].astype(str)
    cov = pcs.set_index(pcs["IID"].astype(str)).loc[order, [f"PC{i}" for i in range(1, cfg["k"] + 1)]].to_numpy()
    return geno, cov


def study_sumstats(panel, trait, source_column):
    ss = read_sumstats(DATA / f"atlas/gwas/{panel}/gwas_{trait}.regenie", sep=" ")
    if panel == "RDP1":
        import importlib.util

        spec = importlib.util.spec_from_file_location("g3", ROOT / "atlas/03_gwas_passports.py")
        g3 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(g3)
        ss = g3.rdp1_published(ss, source_column)
    return ss


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    traits = pd.read_csv(RES / "traits.tsv", sep="\t")
    rows, regions_out = [], []
    for panel, cfg in PANELS.items():
        geno, cov = panel_with_pcs(cfg)
        for t in traits[traits["panel"] == panel].itertuples():
            f = PASS / f"gwas_{panel}_{t.trait}.jsonl"
            if not f.exists():
                continue
            ps = read_passports(f)
            if not ps:
                write_passports(OUT / f.name, ps)
                continue
            t0 = time.time()
            ss = study_sumstats(panel, t.trait, t.source_column)
            n = int(t.n)
            n_signals = 0
            for chrom, a, b, members in loci_regions(ps, cfg["build"]):
                res = finemap_region(ss, chrom, a, b, geno, n=n, L=L, panel_name=cfg["ld"], covariates=cov)
                if res is None:
                    continue
                fm, m, R = res
                n_signals += len(fm.signals)
                sizes = [len(s.variants) for s in fm.signals]
                regions_out.append({"panel": panel, "trait": t.trait, "chrom": chrom, "start": a, "end": b,
                                    "loci": len(members), "snps": len(m), "signals": len(fm.signals),
                                    "cs_sizes": ",".join(map(str, sizes)), "converged": fm.converged,
                                    "ld_outliers": int((np.abs(m["resid"]) > 4).sum())})
                for p in members:
                    p.fine_mapping = to_fine_mapping(fm, m, a, b, L, [cfg["ld"]], lead_pos=p.placement(cfg["build"]).lead_pos, R=R)
            write_passports(OUT / f.name, ps)
            with_signal = sum(p.fine_mapping is not None and p.fine_mapping.lead_signal is not None for p in ps)
            rows.append({"panel": panel, "trait": t.trait, "loci_clumping": len(ps), "signals_finemap": n_signals,
                         "loci_assigned_to_a_signal": with_signal, "seconds": round(time.time() - t0, 1)})
            print(rows[-1], flush=True)
    s = pd.DataFrame(rows)
    s.to_csv(RES / "finemap_loci_summary.tsv", sep="\t", index=False)
    r = pd.DataFrame(regions_out)
    r.to_csv(RES / "finemap_loci_regions.tsv", sep="\t", index=False)
    cs = [int(x) for v in r["cs_sizes"] for x in str(v).split(",") if x]
    print("\n", s[s["loci_clumping"] > 0].to_string(index=False))
    print(f"\nregions {len(r)}; signals {len(cs)}; credible-set size median {np.median(cs):.0f}, "
          f"single-SNP sets {sum(c == 1 for c in cs)}; regions with LD outliers {(r['ld_outliers'] > 0).sum()}")


if __name__ == "__main__":
    main()
