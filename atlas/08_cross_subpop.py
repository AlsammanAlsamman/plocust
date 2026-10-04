"""Atlas step 8: cross-subpopulation fine-mapping on real 3K traits, with dense SNPs.

For traits with strong signals: dense 3K SNPs (4.8M set) around each locus region; REGENIE in the full panel and
separately in subpopulation A (PC1 > 0) and B (PC1 <= 0); each region fine-mapped four ways:
full panel, A alone, B alone, A + B jointly (shared causal variants, own LD, LD PC-adjusted).
Scored on signals, credible-set size of the strongest signal, and whether it overlaps the trait-matched cloned gene.
"""

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from plocust.coloc import to_panel
from plocust.finemap import region_z, susie_rss, susie_rss_multi
from plocust.io import Genotypes, read_sumstats
from plocust.ld import _standardize

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
W = DATA / "atlas" / "xsub"
RES = ROOT / "atlas" / "results"
REGENIE = DATA / "bin/regenie_v4.1.3.gz_x86_64_Linux"
PLINK = DATA / "bin/plink2"
TRAITS = ["culm_length", "glutinous_endosperm", "seed_coat_colour", "apiculus_colour", "hull_colour",
          "leaf_sheath_colour", "auricle_colour", "collar_colour", "lemma_pubescence", "flag_leaf_angle",
          "secondary_branching", "panicle_shattering"]
MAX_SNPS = 3000


def run(cmd, log):
    with open(log, "w") as fh:
        subprocess.run([str(c) for c in cmd], check=True, stdout=fh, stderr=subprocess.STDOUT)


def residual_ld(g, cov):
    x = _standardize(g)
    c = np.column_stack([np.ones(len(x)), cov])
    x = x - c @ np.linalg.lstsq(c, x, rcond=None)[0]
    sd = x.std(0)
    x = x / np.where(sd == 0, 1, sd)
    R = x.T @ x / len(x)
    np.fill_diagonal(R, 1.0)
    return R


def main():
    W.mkdir(parents=True, exist_ok=True)
    traits = pd.read_csv(RES / "traits.tsv", sep="\t").set_index("trait")
    regions = pd.read_csv(RES / "finemap_loci_regions.tsv", sep="\t")
    regions = regions[(regions["panel"] == "3K") & regions["trait"].isin(TRAITS) & (regions["signals"] > 0)]
    known = pd.read_csv(RES / "learn_B_known_genes.tsv", sep="\t")
    print(f"{len(regions)} regions in {regions['trait'].nunique()} traits", flush=True)

    # 1. dense SNPs in all regions
    bed = W / "regions.txt"
    regions[["chrom", "start", "end"]].assign(name="r").to_csv(bed, sep="\t", header=False, index=False)
    dense = W / "dense_loci"
    if not Path(f"{dense}.bed").exists():
        run([PLINK, "--bfile", DATA / "3k/filt/base_filtered_v0.7", "--extract", "bed1", bed, "--maf", "0.01",
             "--make-bed", "--out", dense, "--threads", "8"], W / "extract.log")

    # 2. GWAS: full panel and each subpopulation (REGENIE; step 1 on the pruned core set)
    pheno = DATA / "atlas/gwas/3K/pheno.tsv"
    keep = {"full": None, "A": DATA / "sim/keep_A.txt", "B": DATA / "sim/keep_B.txt"}
    for name, kf in keep.items():
        for binary in (False, True):
            cols = [t for t in TRAITS if bool(traits.loc[t, "binary"]) == binary]
            tag = f"{name}_{'bt' if binary else 'qt'}"
            if all((W / f"{name}_{t}.regenie").exists() for t in cols):
                continue
            k = ["--keep", kf] if kf else []
            bt = ["--bt"] if binary else []
            common = ["--phenoFile", pheno, "--phenoColList", ",".join(cols), "--covarFile", DATA / "3k/covar5.tsv",
                      "--threads", "8", *k, *bt]
            step1 = DATA / "3k/step1set"
            if kf:  # step-1 SNPs must vary inside the subpopulation
                run([PLINK, "--bfile", step1, "--keep", kf, "--maf", "0.05", "--make-bed", "--out", W / f"step1_{name}"],
                    W / f"step1_{name}.log")
                step1 = W / f"step1_{name}"
            try:
                run([REGENIE, "--step", "1", "--bed", step1, "--bsize", "1000", "--lowmem", "--lowmem-prefix",
                     W / f"tmp_{tag}", "--out", W / f"s1_{tag}", *common], W / f"s1_{tag}.log")
                firth = ["--firth", "--approx", "--pThresh", "0.01"] if binary else []
                run([REGENIE, "--step", "2", "--bed", dense, "--pred", W / f"s1_{tag}_pred.list", "--bsize", "400",
                     "--minMAC", "10", "--out", W / name, *firth, *common], W / f"s2_{tag}.log")
            except subprocess.CalledProcessError:
                print(f"REGENIE failed for {tag} (see {W}/s*_{tag}.log)", flush=True)
            print(f"GWAS {tag} done", flush=True)

    # 3. fine-map each region four ways
    geno = Genotypes.open(dense)
    pcs = pd.read_csv(DATA / "3k/pca3k.eigenvec", sep="\t").rename(columns={"#FID": "FID"}).set_index("IID")
    sets = {}
    for name in ("full", "A", "B"):
        ids = pcs.index if name == "full" else pd.read_csv(keep[name], sep="\t", header=None)[1]
        g = geno.subset(ids)
        order = g.samples[g._keep]["iid"]
        sets[name] = (g, pcs.loc[order, [f"PC{i}" for i in range(1, 6)]].to_numpy())
    from plocust.io import read_genes

    g = read_genes(DATA / "ref/Oryza_sativa.IRGSP-1.0.63.gff3.gz").drop_duplicates("id").set_index("id")
    gene_at = {k: {3: int(v["start"]), 4: int(v["end"])} for k, v in g.iterrows()}

    rows = []
    for r in regions.itertuples():
        n = int(traits.loc[r.trait, "n"])
        zs, Rs, ns, matched = {}, {}, {}, {}
        for name, (g, cov) in sets.items():
            f = W / f"{name}_{r.trait}.regenie"
            if not f.exists():
                continue
            ss = read_sumstats(f, sep=" ")
            df, signed = region_z(ss, str(r.chrom), r.start, r.end)
            if not signed or len(df) < 10:
                continue
            m = to_panel(df[["chrom", "pos", "ref", "alt", "z"]].assign(chrom=str(r.chrom)), g, tol=0)
            matched[name] = m.set_index("idx")
            ns[name] = int(ss["n"].median()) if ss["n"].notna().any() else n
        if not {"full", "A", "B"} <= set(matched):
            continue
        common = matched["full"].index.intersection(matched["A"].index).intersection(matched["B"].index)
        if len(common) < 20:
            continue
        # cap: keep the SNPs most significant in any of the three analyses
        score = np.max([np.abs(matched[k].loc[common, "z_a1"].to_numpy()) for k in matched], axis=0)
        common = common[np.sort(np.argsort(-score)[:MAX_SNPS])]
        pos = matched["full"].loc[common, "pos"].to_numpy()
        for name, (g, cov) in sets.items():
            gg = g.read(common.to_numpy())
            Rs[name] = residual_ld(gg, cov)
            zs[name] = matched[name].loc[common, "z_a1"].to_numpy()
        kg = known[(known["panel"] == "3K") & (known["trait"] == r.trait) & (known["chrom"].astype(str) == str(r.chrom))
                   & known["nearest_known_id"].notna() & known["lead"].between(r.start, r.end)]["nearest_known_id"]
        kg = set(kg)
        fits = {
            "full panel": susie_rss(zs["full"], Rs["full"], ns["full"], L=10),
            "A alone": susie_rss(zs["A"], Rs["A"], ns["A"], L=10),
            "B alone": susie_rss(zs["B"], Rs["B"], ns["B"], L=10),
            "A + B joint": susie_rss_multi([zs["A"], zs["B"]], [Rs["A"], Rs["B"]], [ns["A"], ns["B"]], L=10),
        }
        zmax = np.max([np.abs(zs[k]) for k in ("full", "A", "B")], axis=0)
        for label, fm in fits.items():
            sig = [s for s in fm.signals if zmax[s.variants].max() >= 4.42]
            top = max(sig, key=lambda s: s.log_bf) if sig else None
            in_gene = None
            if top is not None and kg:
                cs_pos = pos[top.variants]
                in_gene = any(((cs_pos >= gene_at[k][3] - 5000) & (cs_pos <= gene_at[k][4] + 5000)).any()
                              for k in kg if k in gene_at)
            rows.append({"trait": r.trait, "region": f"{r.chrom}:{r.start}-{r.end}", "snps": len(common),
                         "method": label, "signals": len(sig), "top_cs_size": len(top.variants) if top else None,
                         "known_gene": ";".join(sorted(kg)) or None, "top_cs_hits_known_gene": in_gene})
        print(f"{r.trait} {r.chrom}:{r.start}: {len(common)} SNPs, known={sorted(kg)}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(RES / "xsub_regions.tsv", sep="\t", index=False)
    s = df.groupby("method", sort=False).agg(
        regions=("region", "size"), signals=("signals", "sum"),
        median_top_cs=("top_cs_size", "median"), regions_with_known_gene=("known_gene", lambda x: int(x.notna().sum())),
        top_cs_hits_known_gene=("top_cs_hits_known_gene", lambda x: int((x == True).sum())))  # noqa: E712
    s.to_csv(RES / "xsub_summary.tsv", sep="\t")
    print(s.to_string())


if __name__ == "__main__":
    main()
