"""End-to-end: identify -> anchor -> match -> database -> card, on the synthetic world."""

import numpy as np
import pandas as pd
import pytest

from plocust.anchor import Aligner, anchor_passports, check_uniqueness
from plocust.cli import main
from plocust.db import LocusDB, create, gene_passports
from plocust.identify import IdentifyConfig, StudyInfo, identify_loci
from plocust.io import FlankTable, Genome, Genotypes, normalize_chrom, read_genes, read_sumstats
from plocust.match import MatchConfig, compare, coordinate_overlap
from plocust.passport import AnchorRole, PlacementFlag, read_passports

from .conftest import BLOCK_SNPS, CAUSAL_POS, INSERT_LEN, SNP_STEP

STUDY = StudyInfo(species="Oryza sativa", trait="plant height", study="synthetic", build="SRC")


@pytest.fixture(scope="module")
def loci(world):
    ss = read_sumstats(world.sumstats_path)
    return identify_loci(ss, STUDY, genome=Genome(world.genome, "SRC"), genes=read_genes(world.gff),
                         geno=Genotypes.open(world.bfile), cfg=IdentifyConfig(p_threshold=1e-8))


def test_normalize_chrom():
    assert normalize_chrom("chr01") == normalize_chrom("Chr1") == normalize_chrom(1) == "1"
    assert normalize_chrom("chrUn") == "Un"


def test_read_sumstats_aliases(world):
    ss = read_sumstats(world.sumstats_path)
    assert list(ss.columns[:3]) == ["chrom", "pos", "id"]
    assert ss["chrom"].unique().tolist() == ["1"]
    assert ss["effect_allele"].notna().all()


def test_genotype_reader_roundtrip(world):
    geno = Genotypes.open(world.bfile)
    g = geno.read([0, 1, 2])
    assert g.shape == (300, 3)
    assert set(np.unique(g)) <= {0.0, 2.0}


def test_identify_finds_causal_block(loci):
    assert len(loci) == 1
    p = loci[0]
    blk = p.ld_block
    block_start = (CAUSAL_POS // (SNP_STEP * BLOCK_SNPS)) * SNP_STEP * BLOCK_SNPS + SNP_STEP // 2
    assert blk.start <= CAUSAL_POS <= blk.end
    assert blk.start >= block_start - SNP_STEP and blk.end <= block_start + SNP_STEP * BLOCK_SNPS
    assert 2 <= blk.n_haplotypes <= 4
    assert any(v.pos == CAUSAL_POS for v in p.credible_set.variants)
    assert p.genes[0].id == "G1" and p.genes[0].in_credible_set  # lead may be any SNP in perfect LD with the causal
    roles = {a.role for a in p.anchors}
    assert AnchorRole.lead_snp in roles and len(roles) >= 2  # an edge equal to the lead gets no extra anchor
    lead = p.lead_anchor
    assert lead.sequence[lead.variant_offset] == p.signal.lead.ref
    assert p.notes is None  # lead allele matched the genome


def test_wrong_build_is_flagged(world):
    genome = Genome(world.genome, "SRC")
    ss = read_sumstats(world.sumstats_path)
    # alleles that cannot match the reference base, as when coordinates come from another build
    other = {b: [x for x in "ACGT" if x != b] for b in "ACGT"}
    bases = [genome.fetch("1", int(p), int(p)) for p in ss["pos"]]
    ss["ref"] = [other[b][0] for b in bases]
    ss["alt"] = [other[b][1] for b in bases]
    ps = identify_loci(ss, STUDY, genome=genome, cfg=IdentifyConfig(p_threshold=1e-8))
    assert ps and all("wrong genome build" in (p.notes or "") for p in ps)


def test_anchor_moves_with_insertion(world, loci):
    aligner = Aligner(world.target, build="TGT")
    placed = anchor_passports(loci, aligner)
    pl = placed[0]
    src = loci[0].placement("SRC")
    assert pl.chrom == "1" and pl.strand == "+"
    assert pl.lead_pos == src.lead_pos + INSERT_LEN
    assert pl.anchors_placed == pl.anchors_total == len(loci[0].anchors) and not pl.flags
    assert loci[0].placement("TGT") == pl


def test_anchor_on_inverted_assembly(world, loci):
    aligner = Aligner(world.target_inv, build="INV")
    pl = anchor_passports(loci, aligner)[0]
    src = loci[0].placement("SRC")
    assert pl.strand == "-" and PlacementFlag.inverted in pl.flags
    assert pl.lead_pos == 300_000 - src.lead_pos + 1
    assert PlacementFlag.out_of_order not in pl.flags


def test_unique_anchors(world, loci):
    check_uniqueness(loci, Aligner(world.genome, build="SRC"))
    assert all(a.unique for a in loci[0].anchors)


def test_same_signal_and_distinct(world, loci):
    geno = Genotypes.open(world.bfile)
    a = loci[0]
    # same locus re-detected in a "second study" = same passport data, compared with itself
    pairs = compare([a], [a.model_copy(deep=True)], "SRC", geno)
    assert pairs["call"].tolist() == ["same"]

    # a fake distinct signal in the neighbouring (independent) block
    other = a.model_copy(deep=True)
    other.placements[0].lead_pos += SNP_STEP * BLOCK_SNPS
    other.credible_set = None
    row = compare([a], [other], "SRC", geno, MatchConfig()).iloc[0]
    assert row["call"] == "distinct_nearby"
    assert row["distance_bp"] == SNP_STEP * BLOCK_SNPS
    # the coordinate baseline cannot tell them apart once intervals are padded
    assert coordinate_overlap(a, other, "SRC", pad_kb=50)


def test_gene_database(world, loci, tmp_path):
    genome = Genome(world.genome, "SRC")
    table = pd.DataFrame({"Symbol": ["sd1|GA20ox2", "G2x"], "RAPdb": ["G1", "G2"]})
    kw = pd.DataFrame({"RAPdb": ["G1", "G1"], "Keyword": ["plant height", "gibberellin"]})
    genes = gene_passports(read_genes(world.gff), table, genome, "Oryza sativa", "toy", keywords=kw)
    assert [g.genes[0].name for g in genes] == ["sd1", "G2x"]
    anchor_passports(genes, Aligner(world.target, build="TGT"))
    path = tmp_path / "db.sqlite"
    assert create(path, genes, name="toy", version="0") == 2
    db = LocusDB(path)
    assert db.info()["n_passports"] == "2"
    assert {p.genes[0].name for p in db.by_trait("height")} == {"sd1"}
    near = db.near(loci, "SRC", pad_kb=10)
    pairs = compare(loci, near, "SRC")
    assert pairs.set_index("target")["call"].map(str).to_dict() == {genes[0].passport_id: "gene_in_locus"}
    # the same gene is found on the second assembly
    assert db.get(genes[0].passport_id).placement("TGT").lead_pos > 100_000


def test_flank_table():
    t = FlankTable(pd.DataFrame({"id": ["x"], "seq": ["ACGT[A/G]TTGA"]}), "id", "seq")
    assert t.flank("x") == ("ACGTATTGA", 4)
    assert t.flank("x", "G") == ("ACGTGTTGA", 4)


def test_cli_end_to_end(world, tmp_path):
    out = tmp_path / "loci.jsonl"
    assert main(["identify", "--sumstats", str(world.sumstats_path), "--genome", str(world.genome), "--build", "SRC",
                 "--gff", str(world.gff), "--bfile", str(world.bfile), "--p", "1e-8", "--species", "Oryza sativa",
                 "--trait", "plant height", "--study", "cli", "--check-unique", "--out", str(out)]) == 0
    assert len(read_passports(out)) == 1
    assert (tmp_path / "loci.tsv").exists()
    moved = tmp_path / "moved.jsonl"
    assert main(["anchor", str(out), "--genome", str(world.target), "--build", "TGT", "--out", str(moved)]) == 0
    assert main(["compare", str(out), str(moved), "--build", "TGT", "--out", str(tmp_path / "c.tsv")]) == 0
    assert main(["validate", str(moved)]) == 0
    png = tmp_path / "card.png"
    assert main(["card", str(out), "--sumstats", str(world.sumstats_path), "--gff", str(world.gff), "--bfile",
                 str(world.bfile), "--p", "1e-8", "--out", str(png)]) == 0
    assert png.stat().st_size > 10_000
