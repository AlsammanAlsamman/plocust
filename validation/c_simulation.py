"""Validation C: simulated traits on real 3K genotypes, split into two genetically distinct panels.

Per replicate:
    8 shared causal SNPs           active in both panels           -> loci should be called "same"
    4 nearby-distinct pairs        100-400 kb apart, r2 < 0.1;     -> loci should NOT be called "same"
                                   one SNP active in each panel       (coordinate overlap is fooled here)
    2 private causal SNPs/panel    elsewhere
Each causal SNP explains 4% of phenotypic variance. GWAS per panel: REGENIE mixed model + 5 PCs
(PCs alone leave strong inflation in rice; all replicates are run as one multi-phenotype REGENIE job).

Truth for a detected locus: the active causal SNP (within 1 Mb) with the highest r2 to its lead in its own
panel, if r2 >= 0.3. A pair of loci (panel A x panel B, leads within 1 Mb) is truly "same" when both are
assigned to the same shared causal SNP.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from plocust.identify import IdentifyConfig, StudyInfo, identify_loci
from plocust.io import Genome, Genotypes, read_sumstats
from plocust.ld import r2_matrix
from plocust.match import MatchConfig, compare

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WORK = DATA / "sim"
OUT = ROOT / "validation" / "results" / "c_simulation"
PLINK = DATA / "bin/plink2"
BUILD = "IRGSP-1.0"
N_REP = int(sys.argv[1]) if len(sys.argv) > 1 else 20
H2_EACH = 0.04
REGENIE = DATA / "bin/regenie_v4.1.3.gz_x86_64_Linux"
RNG = np.random.default_rng(2026)


def plink(*args):
    subprocess.run([str(PLINK), "--threads", "8", "--memory", "8000", *map(str, args)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run(cmd):
    subprocess.run([str(x) for x in cmd], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def split_panels():
    pcs = pd.read_csv(DATA / "3k/pca3k.eigenvec", sep="\t").rename(columns={"#FID": "FID"})
    a = pcs[pcs["PC1"] > 0][["FID", "IID"]]
    b = pcs[pcs["PC1"] <= 0][["FID", "IID"]]
    for name, df in (("A", a), ("B", b)):
        df.to_csv(WORK / f"keep_{name}.txt", sep="\t", index=False, header=False)
        if not (WORK / f"freq_{name}.afreq").exists():
            plink("--bfile", DATA / "3k/core3k", "--keep", WORK / f"keep_{name}.txt", "--freq",
                  "--out", WORK / f"freq_{name}")
    return a["IID"].tolist(), b["IID"].tolist(), pcs


def candidates() -> pd.DataFrame:
    fa = pd.read_csv(WORK / "freq_A.afreq", sep="\t", dtype={"ID": str})
    fb = pd.read_csv(WORK / "freq_B.afreq", sep="\t", dtype={"ID": str})
    f = fa[["ID", "ALT_FREQS"]].merge(fb[["ID", "ALT_FREQS"]], on="ID", suffixes=("_A", "_B"))
    maf = lambda x: np.minimum(x, 1 - x)  # noqa: E731
    f = f[(maf(f["ALT_FREQS_A"]) >= 0.1) & (maf(f["ALT_FREQS_B"]) >= 0.1)]
    return f


def draw_causals(geno_all: Genotypes, cand: pd.DataFrame, ga: Genotypes, gb: Genotypes) -> pd.DataFrame:
    v = geno_all.variants.set_index("id")
    pool = v.loc[v.index.intersection(cand["ID"])]
    rows, used = [], []

    def far(chrom, pos):
        return all(c != chrom or abs(p - pos) > 2_000_000 for c, p in used)

    while len([r for r in rows if r["role"] == "shared"]) < 8:
        s = pool.sample(1, random_state=RNG.integers(1e9)).iloc[0]
        if far(s["chrom"], s["pos"]):
            rows.append({"idx": int(s["idx"]), "chrom": s["chrom"], "pos": int(s["pos"]), "role": "shared",
                         "A": True, "B": True})
            used.append((s["chrom"], int(s["pos"])))
    pairs = 0
    while pairs < 4:
        s = pool.sample(1, random_state=RNG.integers(1e9)).iloc[0]
        if not far(s["chrom"], s["pos"]):
            continue
        near = pool[(pool["chrom"] == s["chrom"]) & ((pool["pos"] - s["pos"]).abs().between(100_000, 400_000))]
        if near.empty:
            continue
        t = near.sample(1, random_state=RNG.integers(1e9)).iloc[0]
        idx = [int(s["idx"]), int(t["idx"])]
        if max(r2_matrix(g.read(idx)[:, :1], g.read(idx)[:, 1:])[0, 0] for g in (ga, gb)) >= 0.1:
            continue
        rows.append({"idx": idx[0], "chrom": s["chrom"], "pos": int(s["pos"]), "role": f"pair{pairs}",
                     "A": True, "B": False})
        rows.append({"idx": idx[1], "chrom": t["chrom"], "pos": int(t["pos"]), "role": f"pair{pairs}",
                     "A": False, "B": True})
        used.append((s["chrom"], int(s["pos"])))
        pairs += 1
    for panel in "AB":
        k = 0
        while k < 2:
            s = pool.sample(1, random_state=RNG.integers(1e9)).iloc[0]
            if far(s["chrom"], s["pos"]):
                rows.append({"idx": int(s["idx"]), "chrom": s["chrom"], "pos": int(s["pos"]),
                             "role": f"private{panel}", "A": panel == "A", "B": panel == "B"})
                used.append((s["chrom"], int(s["pos"])))
                k += 1
    return pd.DataFrame(rows)


def phenotype(g: Genotypes, causals: pd.DataFrame, panel: str) -> np.ndarray:
    act = causals[causals[panel]]
    x = g.read(act["idx"].to_numpy())
    x = np.where(np.isnan(x), np.nanmean(x, axis=0), x)
    x = (x - x.mean(0)) / x.std(0)
    y = x @ np.full(len(act), np.sqrt(H2_EACH)) + RNG.normal(0, np.sqrt(1 - H2_EACH * len(act)), len(x))
    return y


def assign(loci, causals, panel, g: Genotypes) -> dict:
    """passport_id -> causal row index (or None)."""
    act = causals[causals[panel]]
    out = {}
    for p in loci:
        pl = p.placement(BUILD)
        near = act[(act["chrom"] == pl.chrom) & ((act["pos"] - pl.lead_pos).abs() <= 1_000_000)]
        lead = g.nearest(pl.chrom, pl.lead_pos, 0)
        best, best_r2 = None, 0.3
        for i, c in near.iterrows():
            r2 = r2_matrix(*(g.read([int(lead["idx"]), int(c["idx"])])[:, k:k + 1] for k in (0, 1)))[0, 0]
            if r2 >= best_r2:
                best, best_r2 = i, r2
        out[p.passport_id] = best
    return out


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    ia, ib, pcs = split_panels()
    geno_all = Genotypes.open(DATA / "3k/core3k")
    ga, gb = geno_all.subset(ia), geno_all.subset(ib)
    genome = Genome(DATA / "ref/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa", BUILD)
    cand = candidates()
    print(f"panel A {ga.n_samples}, panel B {gb.n_samples}, candidate SNPs {len(cand)}", flush=True)

    # 1. draw causal SNPs and phenotypes for every replicate
    reps = [draw_causals(geno_all, cand, ga, gb) for _ in range(N_REP)]
    for panel, g in (("A", ga), ("B", gb)):
        ids = g.samples[g._keep]["iid"].to_numpy()
        ph = pd.DataFrame({"FID": ids, "IID": ids})
        for rep, causals in enumerate(reps):
            ph[f"y{rep}"] = phenotype(g, causals, panel)
        ph.to_csv(WORK / f"pheno_{panel}.tsv", sep="\t", index=False)
    pd.concat({i: c for i, c in enumerate(reps)}).to_csv(OUT / "causals.tsv", sep="\t")

    # 2. one REGENIE run per panel for all replicates
    for panel in "AB":
        if (WORK / f"s2_{panel}_y{N_REP - 1}.regenie").exists() and (WORK / f"pheno_{panel}.tsv").exists():
            print(f"REGENIE panel {panel}: reusing results (replicates are seeded)", flush=True)
            continue
        common = ["--phenoFile", WORK / f"pheno_{panel}.tsv", "--covarFile", DATA / "3k/covar5.tsv",
                  "--keep", WORK / f"keep_{panel}.txt", "--threads", "8"]
        # step-1 SNPs must vary inside the panel
        plink("--bfile", DATA / "3k/step1set", "--keep", WORK / f"keep_{panel}.txt", "--maf", "0.05",
              "--make-bed", "--out", WORK / f"step1_{panel}")
        run([REGENIE, "--step", "1", "--bed", WORK / f"step1_{panel}", "--bsize", "1000", "--lowmem",
             "--lowmem-prefix", WORK / f"tmp{panel}", "--out", WORK / f"s1_{panel}", *common])
        run([REGENIE, "--step", "2", "--bed", DATA / "3k/core3k", "--pred", WORK / f"s1_{panel}_pred.list",
             "--bsize", "400", "--minMAC", "50", "--out", WORK / f"s2_{panel}", *common])
        print(f"REGENIE panel {panel} done", flush=True)

    # 3. passports and matching per replicate
    all_pairs, detect = [], []
    for rep, causals in enumerate(reps):
        loci = {}
        for panel, g in (("A", ga), ("B", gb)):
            ss = read_sumstats(WORK / f"s2_{panel}_y{rep}.regenie", sep=" ")
            loci[panel] = identify_loci(ss, StudyInfo("Oryza sativa", "sim", f"sim{rep}{panel}", BUILD),
                                        genome=genome, geno=g, cfg=IdentifyConfig(window_kb=500, min_significant=2))
        truth_a = assign(loci["A"], causals, "A", ga)
        truth_b = assign(loci["B"], causals, "B", gb)
        for panel, t in (("A", truth_a), ("B", truth_b)):
            hit = {causals.loc[i, "role"] for i in t.values() if i is not None}
            detect.append({"rep": rep, "panel": panel, "loci": len(t),
                           "false_loci": sum(i is None for i in t.values()),
                           "shared_detected": sum(1 for i in set(t.values()) if i is not None
                                                  and causals.loc[i, "role"] == "shared"),
                           "roles": ",".join(sorted(hit))})
        pairs = compare(loci["A"], loci["B"], BUILD, gb, MatchConfig(max_distance_kb=1000, proxy_kb=0))
        if pairs.empty:
            continue
        pa = {p.passport_id: p for p in loci["A"]}
        pb = {p.passport_id: p for p in loci["B"]}
        ca, cb = pairs["query"].map(truth_a), pairs["target"].map(truth_b)
        pairs["truth_same"] = (ca == cb) & ca.notna()
        pairs["truth_role"] = [causals.loc[i, "role"] if i is not None and not pd.isna(i) else "none" for i in ca]
        pairs["role_b"] = [causals.loc[i, "role"] if i is not None and not pd.isna(i) else "none" for i in cb]
        pairs["overlap0"] = [pa[q].placement(BUILD).start <= pb[t].placement(BUILD).end
                             and pb[t].placement(BUILD).start <= pa[q].placement(BUILD).end
                             for q, t in zip(pairs["query"], pairs["target"])]
        pairs["rep"] = rep
        all_pairs.append(pairs)
        print(f"rep {rep}: A {len(loci['A'])} loci, B {len(loci['B'])} loci, pairs {len(pairs)}, "
              f"true same {int(pairs['truth_same'].sum())}", flush=True)

    pairs = pd.concat(all_pairs, ignore_index=True)
    pairs.to_csv(OUT / "pairs.tsv", sep="\t", index=False)
    pd.DataFrame(detect).to_csv(OUT / "detection.tsv", sep="\t", index=False)
    methods = {
        "plocust same-signal (r2>=0.5, panel B LD)": pairs["call"] == "same",
        "coordinate overlap (LD-block intervals)": pairs["overlap0"],
        "lead distance <= 100 kb": pairs["distance_bp"] <= 100_000,
        "lead distance <= 250 kb": pairs["distance_bp"] <= 250_000,
        "lead distance <= 500 kb": pairs["distance_bp"] <= 500_000,
        "colocalization (PP.H4 >= 0.8)": pairs["coloc_call"] == "same",
        "combined: coloc if conclusive, else LD call": pairs["coloc_call"].eq("same")
        | (pairs["coloc_call"].isna() | pairs["coloc_call"].eq("inconclusive")) & pairs["call"].eq("same"),
    }
    truth = pairs["truth_same"]
    trap = pairs["truth_role"].str.startswith("pair") & pairs["role_b"].str.startswith("pair") & ~truth
    rows = []
    for name, pred in methods.items():
        tp = int((pred & truth).sum())
        rows.append({"method": name, "true_same_pairs": int(truth.sum()), "called_same": int(pred.sum()),
                     "precision": round(tp / max(pred.sum(), 1), 3), "recall": round(tp / max(truth.sum(), 1), 3),
                     "nearby_distinct_pairs": int(trap.sum()),
                     "false_same_on_nearby_distinct": int((pred & trap).sum())})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    print(res.to_string(index=False))
    det = pd.DataFrame(detect)
    print(f"\ncoloc calls: {pairs['coloc_call'].value_counts(dropna=False).to_dict()}; "
          f"direction: {pairs['direction'].value_counts(dropna=False).to_dict()}")
    print(f"\nloci per panel-replicate {det['loci'].mean():.1f}; unassigned (false) loci {det['false_loci'].sum()} "
          f"of {det['loci'].sum()}; ambiguous calls {(pairs['call'] == 'ambiguous').sum()}")


if __name__ == "__main__":
    main()
