"""Atlas step 9: do fine-mapped credible sets land on known cloned genes? (RiceVarMap 4K genotypes, incl. deletions)

529 accessions of the RiceVarMap 533 panel with 10 agronomic traits; genotypes from the imputed 4K set
(17.4M SNPs, indels and deletions). For each trait, the regions (gene +- 300 kb) of its classic cloned genes are
GWAS'd (REGENIE) and fine-mapped (SuSiE, PC-adjusted in-sample LD). Scored per gene and trait:
  lead  : is the GWAS lead variant within the gene (+- 5 kb)?
  cs    : does the strongest credible set include a variant within the gene (+- 5 kb)?
"""

import csv
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from plocust.finemap import finemap_region
from plocust.io import Genotypes, read_genes, read_sumstats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
W = DATA / "atlas" / "rvm4k"
RES = ROOT / "atlas" / "results"
REGENIE = DATA / "bin/regenie_v4.1.3.gz_x86_64_Linux"
PLINK = DATA / "bin/plink2"
BED4K = DATA / "rvm/rice4k_geno_add_del"
FLANK, PAD = 300_000, 5_000

GENES = {  # symbol -> (traits, gene ID; GW5 resolved from its MSU locus)
    "GS3": (["Grain_length", "Grain_weight"], "Os03g0407400"),
    "GW5": (["Grain_width", "Grain_weight"], "MSU:LOC_Os05g09520"),
    "GW2": (["Grain_width", "Grain_weight"], "Os02g0244100"),
    "GL7": (["Grain_length"], "Os07g0603300"),
    "qGL3": (["Grain_length"], "Os03g0646900"),
    "GW8": (["Grain_width"], "Os08g0531600"),
    "TGW6": (["Grain_weight"], "Os06g0623700"),
    "GS5": (["Grain_width"], "Os05g0158500"),
    "Hd1": (["Heading_date"], "Os06g0275000"),
    "Ghd7": (["Heading_date", "Plant_height", "Yield"], "Os07g0261200"),
    "DTH8": (["Heading_date"], "Os08g0174500"),
    "Ehd1": (["Heading_date"], "Os10g0463400"),
    "Hd3a": (["Heading_date"], "Os06g0157700"),
    "Ghd7.1": (["Heading_date"], "Os07g0695100"),
    "sd1": (["Plant_height"], "Os01g0883800"),
    "IPA1": (["Num_panicles", "Num_effective_panicles"], "Os08g0509600"),
}


def run(cmd, log):
    with open(log, "w") as fh:
        subprocess.run([str(c) for c in cmd], check=True, stdout=fh, stderr=subprocess.STDOUT)


def gene_coords():
    genes = read_genes(DATA / "ref/Oryza_sativa.IRGSP-1.0.63.gff3.gz").set_index("id")
    msu = pd.read_csv(DATA / "ref/msu7/locus_brief_info.7.0", sep="\t").groupby("locus").agg(
        chrom=("chr", "first"), start=("start", "min"), end=("stop", "max"))
    out = {}
    for sym, (traits, gid) in GENES.items():
        if gid.startswith("MSU:"):
            m = msu.loc[gid[4:]]
            out[sym] = (str(m["chrom"]).replace("Chr", ""), int(m["start"]), int(m["end"]), gid[4:], traits)
        else:
            g = genes.loc[gid]
            out[sym] = (g["chrom"], int(g["start"]), int(g["end"]), gid, traits)
    return out


def main():
    W.mkdir(parents=True, exist_ok=True)
    coords = gene_coords()
    ph = pd.read_csv(DATA / "rvm/phenos.csv")
    keep = W / "keep.txt"
    pd.DataFrame({"f": ph["id_name"], "i": ph["id_name"]}).to_csv(keep, sep="\t", header=False, index=False)
    pheno = W / "pheno.tsv"
    traits = [c for c in ph.columns if c not in ("id", "id_name")]
    ph.rename(columns={"id_name": "IID"}).assign(FID=lambda d: d["IID"])[["FID", "IID", *traits]].to_csv(
        pheno, sep="\t", index=False, na_rep="NA")

    # genotypes: gene regions (all variants) and a pruned genome-wide set for REGENIE step 1 and PCs
    regions = W / "regions.txt"
    pd.DataFrame([(c, max(1, a - FLANK), b + FLANK, s) for s, (c, a, b, _, _) in coords.items()]).to_csv(
        regions, sep="\t", header=False, index=False)
    if not (W / "regions.bed").exists():
        run([PLINK, "--bfile", BED4K, "--keep", keep, "--extract", "bed1", regions, "--maf", "0.01",
             "--make-bed", "--out", W / "regions", "--threads", "8", "--memory", "6000"], W / "regions.log")
    if not (W / "step1.bed").exists():
        run([PLINK, "--bfile", BED4K, "--keep", keep, "--snps-only", "just-acgt", "--maf", "0.05", "--geno", "0.1",
             "--thin-count", "400000", "--seed", "1", "--make-bed", "--out", W / "thin", "--threads", "8",
             "--memory", "6000"], W / "thin.log")
        run([PLINK, "--bfile", W / "thin", "--indep-pairwise", "200kb", "0.5", "--out", W / "prune"], W / "prune.log")
        run([PLINK, "--bfile", W / "thin", "--extract", W / "prune.prune.in", "--make-bed", "--out", W / "step1"],
            W / "step1.log")
        run([PLINK, "--bfile", W / "step1", "--pca", "5", "--out", W / "pca"], W / "pca.log")
    covar = W / "covar.tsv"
    pd.read_csv(W / "pca.eigenvec", sep="\t").rename(columns={"#FID": "FID"}).to_csv(covar, sep="\t", index=False)
    if not (W / f"gwas_{traits[-1]}.regenie").exists():
        common = ["--phenoFile", pheno, "--phenoColList", ",".join(traits), "--covarFile", covar, "--threads", "8"]
        run([REGENIE, "--step", "1", "--bed", W / "step1", "--bsize", "1000", "--lowmem", "--lowmem-prefix",
             W / "tmp", "--out", W / "s1", *common], W / "s1.log")
        run([REGENIE, "--step", "2", "--bed", W / "regions", "--pred", W / "s1_pred.list", "--bsize", "400",
             "--minMAC", "10", "--out", W / "gwas", *common], W / "s2.log")
    print("GWAS done", flush=True)

    geno = Genotypes.open(W / "regions")
    pcs = pd.read_csv(covar, sep="\t").set_index("IID")
    order = geno.samples["iid"]
    cov = pcs.loc[order, [f"PC{i}" for i in range(1, 6)]].to_numpy()
    rows = []
    for sym, (chrom, a, b, gid, gtraits) in coords.items():
        for trait in gtraits:
            ss = read_sumstats(W / f"gwas_{trait}.regenie", sep=" ")
            reg = ss[(ss["chrom"] == chrom) & ss["pos"].between(a - FLANK, b + FLANK)]
            if reg.empty:
                continue
            lead = reg.loc[reg["p"].idxmin()]
            row = {"gene": sym, "gene_id": gid, "trait": trait, "chrom": chrom, "gene_start": a, "gene_end": b,
                   "min_p": float(lead["p"]), "lead_pos": int(lead["pos"]),
                   "lead_in_gene": bool(a - PAD <= lead["pos"] <= b + PAD),
                   "lead_kb_from_gene": round(max(0, a - lead["pos"], lead["pos"] - b) / 1000, 1)}
            for ld_label, c in (("raw", None), ("pc_adjusted", cov)):
                res = finemap_region(ss, chrom, a - FLANK, b + FLANK, geno, n=529, L=10, covariates=c, max_snps=3000)
                if res is None:
                    continue
                fm, m, _ = res
                in_gene = (m["pos"] >= a - PAD) & (m["pos"] <= b + PAD)
                row[f"{ld_label}_max_pip_in_gene"] = round(float(fm.pip[in_gene.to_numpy()].max()), 3) if in_gene.any() else None
                row[f"{ld_label}_signals"] = len(fm.signals)
                row[f"{ld_label}_cs_in_gene"] = any(bool(in_gene.iloc[s.variants].any()) for s in fm.signals)
                row[f"{ld_label}_ld_outliers"] = int((np.abs(m["resid"]) > 4).sum())
            rows.append(row)
            print({k: row.get(k) for k in ("gene", "trait", "min_p", "lead_in_gene", "raw_max_pip_in_gene",
                                           "pc_adjusted_max_pip_in_gene")}, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "known_genes_4k.tsv", sep="\t", index=False)
    sig = df[df["min_p"] < 1e-5]
    print(f"\n{len(df)} gene x trait tests; {len(sig)} with a signal (p < 1e-5) in the gene region")
    print(f"  lead variant within gene +-5 kb          : {int(sig['lead_in_gene'].sum())}/{len(sig)}")
    for ld in ("raw", "pc_adjusted"):
        print(f"  {ld:11} LD: a credible set hits the gene : {int(sig[f'{ld}_cs_in_gene'].fillna(False).sum())}/{len(sig)}; "
              f"gene variant with PIP >= 0.5: {int((sig[f'{ld}_max_pip_in_gene'] >= 0.5).sum())}/{len(sig)}")


if __name__ == "__main__":
    main()
