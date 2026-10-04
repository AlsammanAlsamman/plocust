"""Validation F: the second round of features on real rice data.

  F1  genome-build detection         RDP1 and 3K summary statistics vs MSU6 / IRGSP-1.0 / MH63
  F2  marker lookup (no genome)      RDP1 passports from chip-marker flanks stored in the database
  F3  cross-study, all evidence      RDP1 (REGENIE, MSU6) vs 3K (REGENIE, IRGSP-1.0) at sd1, Wx, Rc:
                                     LD call (core and dense panels), effect direction, imprint profile,
                                     colocalization
  F4  locus cards with the regional imprint

RDP1 is re-analysed with REGENIE here because the published p-values carry no effect signs,
which the direction test and the imputation step of colocalization need.
"""

import shutil
from pathlib import Path

import pandas as pd

from plocust.anchor import Aligner, anchor_passports
from plocust.build import detect_build
from plocust.db import LocusDB, add_markers, attach_ld_panel
from plocust.identify import IdentifyConfig, StudyInfo, identify_loci
from plocust.io import Genome, Genotypes, read_genes, read_sumstats
from plocust.match import MatchConfig, compare
from plocust.passport import write_passports
from plocust.plot import locus_card

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "validation" / "results" / "f_new_features"
OUT.mkdir(parents=True, exist_ok=True)
IRGSP = "IRGSP-1.0"
SPECIES = "Oryza sativa"

# gene -> (region on IRGSP-1.0, RDP1 trait column, 3K trait, trait name)
GENES = {
    "sd1": (("1", 37_000_000, 39_600_000), "plant_height", "culm_length", "Plant height"),
    "Wx": (("6", 1_000_000, 2_600_000), "amylose", "endosperm", "Amylose content"),
    "Rc": (("7", 5_400_000, 6_800_000), "pericarp", "seed_coat", "Pericarp color"),
}
# MSU6 -> IRGSP-1.0 shifts are small here, so the same windows (+-1 Mb slack) select the RDP1 loci on MSU6
SLACK = 1_000_000


def genomes():
    return {"MSU6": Genome(DATA / "msu6/MSU6.fa", "MSU6"),
            IRGSP: Genome(DATA / "ref/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa", IRGSP),
            "MH63RS2": Genome(DATA / "asm/Oryza_sativa_mh63.MH63RS2.dna.toplevel.fa", "MH63RS2")}


def f1_build_detection(gs):
    rows = []
    for name, path in (("RDP1 (REGENIE, MSU6)", DATA / "rdp1/rs2_plant_height.regenie"),
                       ("3K (REGENIE, IRGSP-1.0)", DATA / "3k/s2_culm_length.regenie")):
        res = detect_build(read_sumstats(path, sep=" "), list(gs.values()))
        res.insert(0, "study", name)
        rows.append(res)
    df = pd.concat(rows)
    df.to_csv(OUT / "f1_build_detection.tsv", sep="\t", index=False)
    return df


def f2_marker_lookup(gs):
    """Markers from the RDP1 chip manifest (20+20 bp MSU7 flanks, and 16+16 bp MSU6 flanks)."""
    db = OUT / "markers.sqlite"
    shutil.copy(DATA / "db/plocust-db-rice.sqlite", db)
    f7 = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU7.SNP_flanking_seq.txt", sep="\t", dtype=str)
    f7["sequence"] = f7["msu7_seq1"] + "[" + f7["alleles"] + "]" + f7["msu7_seq2"]
    f6 = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU6.SNP_flanking_seq.txt", sep="\t", dtype=str)
    f6["sequence"] = f6["X5p_MSU6"] + "[" + f6["alleles"] + "]" + f6["X3p_MSU6"]
    add_markers(db, f7.rename(columns={"snp_id": "marker_id"}), "RDP1 44K (41 bp)")
    info = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU6.SNP_Information.MSU7.txt", sep="\t", index_col=False)
    truth = dict(zip(info["SNPID"], pd.to_numeric(info["position.IRGSP1-MSU7"], errors="coerce")))
    ss = read_sumstats(DATA / "rdp1/rs2_amylose.regenie", sep=" ")
    aligner = Aligner(gs[IRGSP].path, build=IRGSP)
    rows = []
    for label, table in (("41 bp (MSU7 manifest)", f7), ("33 bp (MSU6 manifest)", f6)):
        tmp = OUT / "m.sqlite"
        shutil.copy(DATA / "db/plocust-db-rice.sqlite", tmp)
        add_markers(tmp, table.rename(columns={"snp_id": "marker_id"}), label)
        flanks = LocusDB(tmp).marker_flanks()
        loci = identify_loci(ss, StudyInfo(SPECIES, "Amylose content", "RDP1", "MSU6"), flanks=flanks,
                             cfg=IdentifyConfig(p_threshold=1e-4, min_significant=1))
        anchor_passports(loci, aligner)
        for p in loci:
            pl = p.placement(IRGSP)
            rows.append({"flank": label, "lead": p.signal.lead.id, "anchor_bp": len(p.lead_anchor.sequence),
                         "placed": pl is not None, "exact": pl is not None and pl.lead_pos == truth.get(p.signal.lead.id)})
        tmp.unlink()
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "f2_marker_lookup.tsv", sep="\t", index=False)
    return df.groupby("flank")[["placed", "exact"]].agg(["sum", "size"])


def in_region(p, build, region, slack=0):
    pl = p.placement(build)
    chrom, lo, hi = region
    return pl is not None and pl.chrom == chrom and lo - slack <= pl.lead_pos <= hi + slack


def f3_cross_study(gs):
    genes = read_genes(DATA / "ref/Oryza_sativa.IRGSP-1.0.63.gff3.gz")
    aligner = Aligner(gs[IRGSP].path, build=IRGSP)
    core, dense = Genotypes.open(DATA / "3k/core3k"), Genotypes.open(DATA / "3k/dense_regions")
    g_rdp1 = Genotypes.open(DATA / "rdp1/rdp1")
    db = OUT / "markers.sqlite"
    attach_ld_panel(db, DATA / "3k/dense_regions", IRGSP, "3K filtered v0.7, sd1/Wx/Rc regions")
    rows, cards = [], []
    for gene, (region, rdp_col, k3_col, trait) in GENES.items():
        ss_r = read_sumstats(DATA / f"rdp1/rs2_{rdp_col}.regenie", sep=" ")
        ss_k = read_sumstats(DATA / f"3k/s2_{k3_col}.regenie", sep=" ")
        rdp = identify_loci(ss_r, StudyInfo(SPECIES, trait, "RDP1 (REGENIE)", "MSU6", n_samples=int(ss_r["n"].median())),
                            genome=gs["MSU6"], geno=g_rdp1, cfg=IdentifyConfig(p_threshold=1e-5, min_significant=1))
        rdp = [p for p in rdp if in_region(p, "MSU6", region, SLACK)]
        anchor_passports(rdp, aligner)
        rdp = [p for p in rdp if in_region(p, IRGSP, region)]
        k3 = identify_loci(ss_k, StudyInfo(SPECIES, trait, "3K RG (REGENIE)", IRGSP, n_samples=int(ss_k["n"].median())),
                           genome=gs[IRGSP], genes=genes, geno=core, cfg=IdentifyConfig(min_significant=2))
        k3 = [p for p in k3 if in_region(p, IRGSP, region)]
        write_passports(OUT / f"rdp1_{gene}.jsonl", rdp)
        write_passports(OUT / f"3k_{gene}.jsonl", k3)
        for panel_name, panel in (("core 365K", core), ("dense 4.8M", dense)):
            pairs = compare(rdp, k3, IRGSP, panel, MatchConfig(max_distance_kb=1000, proxy_kb=10))
            pairs.insert(0, "panel", panel_name)
            pairs.insert(0, "gene", gene)
            rows.append(pairs)
        cards.append((gene, k3_col, k3, ss_k))
    df = pd.concat(rows, ignore_index=True)
    df.to_csv(OUT / "f3_pairs.tsv", sep="\t", index=False)

    known_db = LocusDB(DATA / "db/plocust-db-rice.sqlite")
    for gene, k3_col, k3, ss_k in cards:  # F4
        if not k3:
            continue
        best = min(k3, key=lambda p: p.signal.pvalue)
        known = [k for k in known_db.near([best], IRGSP, 150) if k.kind.value == "gene"]
        locus_card(best, OUT / f"card_{gene}.png", sumstats=ss_k, genes=genes, geno=core, known=known,
                   p_threshold=0.05 / len(ss_k))
    return df


if __name__ == "__main__":
    gs = genomes()
    print("F1 BUILD DETECTION\n", f1_build_detection(gs).to_string(index=False), flush=True)
    print("\nF2 MARKER LOOKUP (no genome; anchored onto IRGSP-1.0)\n", f2_marker_lookup(gs).to_string(), flush=True)
    df = f3_cross_study(gs)
    cols = ["gene", "panel", "query_lead", "target_lead", "distance_bp", "r2", "direction", "profile_corr", "call",
            "PP.H3", "PP.H4", "coloc_snps", "coloc_call"]
    best = df.sort_values("PP.H4", ascending=False, na_position="last").groupby(["gene", "panel"]).head(3)
    print("\nF3 CROSS-STUDY (top pairs per gene and panel)\n", best[cols].to_string(index=False))
