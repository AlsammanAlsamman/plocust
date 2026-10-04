"""Validation D + E: the same traits mapped in two real rice panels on two genome builds.

    RDP1  413 accessions, 44K array, published mixed-model p-values, MSU6 coordinates (Zhao et al. 2011)
    3K    ~2,100 phenotyped accessions, 365K core SNPs, our REGENIE mixed-model GWAS, IRGSP-1.0 coordinates

Steps: passports for both studies -> anchor RDP1 from MSU6 onto IRGSP-1.0 -> check the anchored
positions against the published MSU7 (= IRGSP-1.0) positions -> compare loci between studies with
the same-signal test in the 3K panel -> match both studies to the known-gene database.
"""

import sys
from pathlib import Path

import pandas as pd

from plocust.anchor import Aligner, anchor_passports, check_uniqueness
from plocust.db import LocusDB
from plocust.identify import IdentifyConfig, StudyInfo, identify_loci
from plocust.io import Genome, Genotypes, read_genes, read_sumstats, standardize_sumstats
from plocust.match import MatchConfig, compare
from plocust.passport import read_passports, write_passports
from plocust.plot import locus_card

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "validation" / "results" / "d_real_gwas"
OUT.mkdir(parents=True, exist_ok=True)
IRGSP = "IRGSP-1.0"

# RDP1 trait -> (3K trait, Planteome TO term, known genes to look for)
TRAITS = {
    "Plant height": ("culm_length", "TO:0000207", ["sd1"]),
    "Awn presence": ("awn", "TO:0000141", ["An-1", "LABA1", "GAD1", "RAE2", "An-2", "DL"]),
    "Pericarp color": ("seed_coat", "TO:0000698", ["Rc", "Rd"]),
    "Amylose content": ("endosperm", "TO:0000196", ["Wx"]),
}
# Panicle length and flag leaf length were dropped: no RDP1 SNP reaches p < 1e-5 for either.
SPECIES = "Oryza sativa"
CFG_RDP1 = IdentifyConfig(p_threshold=1e-5, window_kb=500, min_significant=1)  # 44K array, 413 lines: low power
CFG_3K = IdentifyConfig(window_kb=500, min_significant=2)  # Bonferroni


def rdp1_sumstats(trait: str) -> pd.DataFrame:
    pv = pd.read_csv(DATA / "rdp1/MixedModel_Pval_all.txt", sep="\t", usecols=["SNPID", "MAF", trait])
    bim = pd.read_csv(DATA / "rdp1/rdp1.bim", sep="\t", header=None, names=["chrom", "id", "cm", "pos", "a1", "a2"])
    m = pv.merge(bim, left_on="SNPID", right_on="id")
    raw = pd.DataFrame({"chrom": m["chrom"], "pos": m["pos"], "id": m["id"], "ref": m["a2"], "alt": m["a1"],
                        "p": m[trait], "af": m["MAF"], "n": 413})
    return standardize_sumstats(raw)


def build_passports(force: bool = False):
    msu6 = Genome(DATA / "msu6/MSU6.fa", "MSU6")
    irgsp = Genome(DATA / "ref/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa", IRGSP)
    genes = read_genes(DATA / "ref/Oryza_sativa.IRGSP-1.0.63.gff3.gz")
    g_rdp1 = Genotypes.open(DATA / "rdp1/rdp1")
    g_3k = Genotypes.open(DATA / "3k/core3k")
    a_irgsp = Aligner(irgsp.path, build=IRGSP)
    a_msu6 = Aligner(msu6.path, build="MSU6")
    for rdp_trait, (k3_trait, to, _) in TRAITS.items():
        tag = k3_trait
        f_rdp, f_3k = OUT / f"rdp1_{tag}.jsonl", OUT / f"3k_{tag}.jsonl"
        if f_rdp.exists() and f_3k.exists() and not force:
            continue
        ss = rdp1_sumstats(rdp_trait)
        rdp = identify_loci(ss, StudyInfo(SPECIES, rdp_trait, "RDP1 (Zhao et al. 2011)", "MSU6", to, "RDP1 413",
                                          413, "EMMA mixed model", "10.1038/ncomms1467"),
                            genome=msu6, geno=g_rdp1, cfg=CFG_RDP1)
        check_uniqueness(rdp, a_msu6)
        anchor_passports(rdp, a_irgsp)
        write_passports(f_rdp, rdp)

        ss3 = read_sumstats(DATA / f"3k/s2_{k3_trait}.regenie", sep=" ")
        k3 = identify_loci(ss3, StudyInfo(SPECIES, rdp_trait, "3K RG (this study)", IRGSP, to, "3K RG core SNPs",
                                          None, "REGENIE mixed model + 5 PCs"),
                           genome=irgsp, genes=genes, geno=g_3k, cfg=CFG_3K)
        check_uniqueness(k3, a_irgsp)
        write_passports(f_3k, k3)
        print(f"{rdp_trait}: RDP1 {len(rdp)} loci, 3K {len(k3)} loci", flush=True)


def anchoring_accuracy() -> pd.DataFrame:
    """Anchored IRGSP-1.0 lead positions vs the published MSU7 positions of the same SNPs."""
    info = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU6.SNP_Information.MSU7.txt", sep="\t", index_col=False)
    msu7 = dict(zip(info["SNPID"], pd.to_numeric(info["position.IRGSP1-MSU7"], errors="coerce")))
    rows = []
    for f in sorted(OUT.glob("rdp1_*.jsonl")):
        for p in read_passports(f):
            src, pl = p.placement("MSU6"), p.placement(IRGSP)
            truth = msu7.get(p.signal.lead.id)
            rows.append({"trait": p.trait.name, "lead": p.signal.lead.id, "chrom": src.chrom, "msu6_pos": src.lead_pos,
                         "msu7_truth": truth, "anchored_pos": pl.lead_pos if pl else None,
                         "anchored_chrom": pl.chrom if pl else None,
                         "flags": ",".join(f.value for f in pl.flags) if pl else "unplaced",
                         "unique_anchor": p.lead_anchor.unique})
    df = pd.DataFrame(rows)
    df["exact"] = (df["anchored_pos"] == df["msu7_truth"]) & (df["anchored_chrom"] == df["chrom"])
    df["build_shift_bp"] = (df["msu7_truth"] - df["msu6_pos"]).abs()
    df.to_csv(OUT / "anchoring_rdp1_msu6_to_irgsp.tsv", sep="\t", index=False)
    return df


def cross_study() -> pd.DataFrame:
    g_3k = Genotypes.open(DATA / "3k/core3k")
    rows, all_pairs = [], []
    for rdp_trait, (tag, _, _) in TRAITS.items():
        rdp, k3 = read_passports(OUT / f"rdp1_{tag}.jsonl"), read_passports(OUT / f"3k_{tag}.jsonl")
        pairs = compare(rdp, k3, IRGSP, g_3k, MatchConfig(max_distance_kb=1000, proxy_kb=10))
        pairs.insert(0, "trait", rdp_trait)
        all_pairs.append(pairs)
        best = pairs.sort_values("score", ascending=False).groupby("query").head(1)
        # naive baseline: MSU6 coordinates used as if they were IRGSP-1.0, lead within 100 kb
        k3_leads = [(p.placement(IRGSP).chrom, p.placement(IRGSP).lead_pos) for p in k3]
        naive = sum(any(c == p.placement("MSU6").chrom and abs(x - p.placement("MSU6").lead_pos) <= 100_000
                        for c, x in k3_leads) for p in rdp)
        anchored_pos = sum(any(c == p.placement(IRGSP).chrom and abs(x - p.placement(IRGSP).lead_pos) <= 100_000
                               for c, x in k3_leads) for p in rdp if p.placement(IRGSP))
        rows.append({"trait": rdp_trait, "rdp1_loci": len(rdp), "3k_loci": len(k3),
                     "rdp1_anchored": sum(p.placement(IRGSP) is not None for p in rdp),
                     "same_signal": int((best["call"] == "same").sum()),
                     "ambiguous": int((best["call"] == "ambiguous").sum()),
                     "distinct_nearby_only": int((best["call"] == "distinct_nearby").sum()),
                     "within100kb_anchored": anchored_pos, "within100kb_naive_msu6": naive})
    pd.concat(all_pairs).to_csv(OUT / "rdp1_vs_3k_pairs.tsv", sep="\t", index=False)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "rdp1_vs_3k_summary.tsv", sep="\t", index=False)
    return summary


def known_genes() -> pd.DataFrame:
    db = LocusDB(DATA / "db/plocust-db-rice.sqlite")
    rows = []
    for rdp_trait, (tag, _, expected) in TRAITS.items():
        for study, f in (("RDP1", OUT / f"rdp1_{tag}.jsonl"), ("3K", OUT / f"3k_{tag}.jsonl")):
            loci = [p for p in read_passports(f) if p.placement(IRGSP)]
            known = db.near(loci, IRGSP, 500)
            pairs = compare(loci, known, IRGSP, cfg=MatchConfig(max_distance_kb=500, gene_pad_kb=0))
            sym = {k.passport_id: k.genes[0].name for k in known}
            pairs["gene"] = pairs["target"].map(sym)
            inside = pairs[pairs["call"] == "gene_in_locus"]
            for gene in expected:
                hit = pairs[pairs["gene"].str.split("|").str[0].str.lower() == gene.lower()] if len(pairs) else pairs
                rows.append({"trait": rdp_trait, "study": study, "expected_gene": gene,
                             "in_locus": bool(len(hit) and (hit["call"] == "gene_in_locus").any()),
                             "nearest_kb": round(hit["distance_bp"].min() / 1000, 1) if len(hit) else None})
            rows.append({"trait": rdp_trait, "study": study, "expected_gene": "(all cloned genes in loci)",
                         "in_locus": len(inside), "nearest_kb": None})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "known_genes.tsv", sep="\t", index=False)
    return df


def cards():
    db = LocusDB(DATA / "db/plocust-db-rice.sqlite")
    genes = read_genes(DATA / "ref/Oryza_sativa.IRGSP-1.0.63.gff3.gz")
    g = Genotypes.open(DATA / "3k/core3k")
    for tag, gene in (("seed_coat", "Rc"), ("endosperm", "Wx"), ("culm_length", "sd1")):
        loci = read_passports(OUT / f"3k_{tag}.jsonl")
        ss = read_sumstats(DATA / f"3k/s2_{tag}.regenie", sep=" ")
        target = [k for k in db.near(loci, IRGSP, 300) if k.genes and k.genes[0].name.lower() == gene.lower()]
        if not target:
            continue
        tl = target[0].placement(IRGSP)
        near = [p for p in loci if p.placement(IRGSP).chrom == tl.chrom
                and abs(p.placement(IRGSP).lead_pos - tl.lead_pos) <= 300_000]
        best = min(near, key=lambda p: p.signal.pvalue)  # strongest locus near the gene
        known = [k for k in db.near([best], IRGSP, 150) if k.kind.value == "gene"]
        locus_card(best, OUT / f"card_3k_{tag}_{gene}.png", sumstats=ss, genes=genes, geno=g, known=known,
                   p_threshold=0.05 / len(ss))


if __name__ == "__main__":
    steps = sys.argv[1:] or ["build", "anchoring", "compare", "genes", "cards"]
    if "build" in steps:
        build_passports()
    if "anchoring" in steps:
        a = anchoring_accuracy()
        placed = a["anchored_pos"].notna()
        print(f"\nANCHORING RDP1 MSU6 -> IRGSP-1.0: {len(a)} lead SNPs, placed {placed.sum()}, "
              f"exact position {a['exact'].sum()} ({a['exact'].mean():.1%}); "
              f"median MSU6->MSU7 shift {a['build_shift_bp'].median():,.0f} bp")
    if "compare" in steps:
        print("\nCROSS-STUDY (best match per RDP1 locus, LD in 3K panel)\n", cross_study().to_string(index=False))
    if "genes" in steps:
        print("\nKNOWN GENES\n", known_genes().to_string(index=False))
    if "cards" in steps:
        cards()
