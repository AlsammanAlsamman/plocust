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


def test_read_sumstats_log10p(tmp_path):
    f = tmp_path / "r.regenie"
    f.write_text("CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N BETA SE LOG10P\n1 100 a C T 0.2 500 0.1 0.02 5\n1 200 b G A 0.3 500 0.1 0.02 400\n")
    ss = read_sumstats(f, sep=" ")
    assert ss["p"].tolist() == pytest.approx([1e-5, 1e-300], rel=1e-9)
    assert ss["pos"].tolist() == [100, 200] and ss["af"].tolist() == [0.2, 0.3]



def test_credible_sets_do_not_borrow_other_clumps(world, tmp_path):
    """A strong signal with a perfect-LD copy at 140 kb, and an independent weaker signal at 125 kb.

    LD clumping gives the 140 kb copy to the strong clump. It lies closer to the weak lead, and
    before the fix it ended up in the weak locus's credible set.
    """
    from .conftest import write_bed

    rng = np.random.default_rng(3)
    n = 400
    x = rng.integers(0, 2, n) * 2.0
    z = rng.integers(0, 2, n) * 2.0
    noise = rng.integers(0, 2, (n, 6)) * 2.0
    g = np.column_stack([x, z, x, noise])
    pos = [100_500, 125_500, 140_500, 20_500, 40_500, 60_500, 200_500, 230_500, 260_500]
    order = np.argsort(pos)
    g, pos = g[:, order], np.array(pos)[order]
    bim = pd.DataFrame({"chrom": "1", "id": [f"v{p}" for p in pos], "cm": 0, "pos": pos, "a1": "N", "a2": "N"})
    write_bed(tmp_path / "p", g, bim, [f"i{k}" for k in range(n)])
    geno = Genotypes.open(tmp_path / "p")
    p_of = {100_500: 1e-40, 140_500: 1e-40 * 1.0001, 125_500: 1e-12}
    ss = pd.DataFrame({"chrom": "1", "pos": pos, "id": [f"v{p}" for p in pos],
                       "p": [p_of.get(p, 0.5) for p in pos]})
    from plocust.io import standardize_sumstats

    ps = identify_loci(standardize_sumstats(ss), STUDY, genome=Genome(world.genome, "SRC"), geno=geno,
                       cfg=IdentifyConfig(p_threshold=1e-8))
    weak = next(p for p in ps if p.signal.lead.pos == 125_500)
    assert 140_500 not in {v.pos for v in weak.credible_set.variants}


def test_detect_build(world):
    from plocust.build import detect_build

    ss = read_sumstats(world.sumstats_path)
    res = detect_build(ss, [Genome(world.genome, "SRC"), Genome(world.target, "TGT")])
    assert res.loc[0, "build"] == "SRC" and res.loc[0, "verdict"] == "likely build"
    assert res.loc[0, "allele_match"] == 1.0
    assert res.loc[1, "score"] < 0.8  # 5 kb insertion shifts every SNP after 50 kb


def test_cli_detect_build(world, capsys):
    assert main(["detect-build", "--sumstats", str(world.sumstats_path), "--genome", f"{world.genome}:SRC",
                 "--genome", f"{world.target}:TGT"]) == 0
    assert "likely build" in capsys.readouterr().out


def test_effect_direction(world, loci):
    geno = Genotypes.open(world.bfile)
    a = loci[0]
    same = a.model_copy(deep=True)
    assert compare([a], [same], "SRC", geno).iloc[0]["direction"] == "concordant"
    flipped = a.model_copy(deep=True)
    flipped.signal.beta = -flipped.signal.beta  # same SNP, opposite effect: cannot be the same causal effect
    row = compare([a], [flipped], "SRC", geno).iloc[0]
    assert row["direction"] == "discordant" and row["call"] == "same_opposite_effect"
    swapped = flipped.model_copy(deep=True)  # opposite sign but reported for the other allele: concordant again
    swapped.signal.effect_allele = swapped.signal.lead.ref
    assert compare([a], [swapped], "SRC", geno).iloc[0]["direction"] == "concordant"


def test_imprint(loci):
    imp = loci[0].imprint
    assert imp.z_signed and len(imp.z) == len(imp.pos) == len(imp.ids) > 100
    assert len(imp.profile) == imp.bins and max(imp.profile) == pytest.approx(max(abs(z) for z in imp.z))
    assert imp.n_significant >= 1 and imp.independent_signals == 1
    assert imp.ld_edges and all(r2 >= 0.5 for _, _, r2 in imp.ld_edges)
    lead = loci[0].signal.lead.pos
    assert imp.pos.index(lead) == int(np.argmax(np.abs(imp.z)))


def test_profile_similarity(world, loci):
    a = loci[0]
    row = compare([a], [a.model_copy(deep=True)], "SRC").iloc[0]
    assert row["profile_corr"] == pytest.approx(1.0)


def test_marker_lookup_without_genome(world, tmp_path):
    """A breeder with marker IDs only: flanks come from the database, then anchors place the locus."""
    from plocust.db import add_markers

    genome = Genome(world.genome, "SRC")
    ss = read_sumstats(world.sumstats_path)
    seqs = []
    for pos, ref, alt in zip(ss["pos"], ss["ref"], ss["alt"]):
        left, right = genome.fetch("1", pos - 60, pos - 1), genome.fetch("1", pos + 1, pos + 60)
        seqs.append(f"{left}[{ref}/{alt}]{right}")
    db = tmp_path / "db.sqlite"
    create(db, [], name="t", version="0")
    assert add_markers(db, pd.DataFrame({"marker_id": ss["id"], "sequence": seqs}), "toy chip") == len(ss)
    flanks = LocusDB(db).marker_flanks()
    ps = identify_loci(ss, STUDY, flanks=flanks, cfg=IdentifyConfig(p_threshold=1e-8))
    assert ps and len(ps[0].lead_anchor.sequence) == 121
    pl = anchor_passports(ps, Aligner(world.target, build="TGT"))[0]
    assert pl.lead_pos == ps[0].signal.lead.pos + INSERT_LEN

    # chip manifests often give only ~20 bp each side: short anchors use the near-perfect unique-hit rule
    short = pd.DataFrame({"marker_id": ss["id"], "sequence": [q[40:-40] for q in seqs]})
    db2 = tmp_path / "db2.sqlite"
    create(db2, [], name="t", version="0")
    add_markers(db2, short, "short chip")
    ps = identify_loci(ss, STUDY, flanks=LocusDB(db2).marker_flanks(), cfg=IdentifyConfig(p_threshold=1e-8))
    assert len(ps[0].lead_anchor.sequence) == 41
    pl = anchor_passports(ps, Aligner(world.target, build="TGT"))[0]
    assert pl.lead_pos == ps[0].signal.lead.pos + INSERT_LEN


def test_reference_ld_panel(world, loci, tmp_path):
    from plocust.db import attach_ld_panel
    from plocust.passport import write_passports

    db = tmp_path / "db.sqlite"
    create(db, [], name="t", version="0")
    attach_ld_panel(db, world.bfile, "SRC", "toy panel")
    g = LocusDB(db).ld_panel("SRC")
    assert g is not None and g.n_samples == 300 and LocusDB(db).ld_panel("OTHER") is None
    f = tmp_path / "l.jsonl"
    write_passports(f, loci)
    out = tmp_path / "c.tsv"
    assert main(["compare", str(f), str(f), "--build", "SRC", "--db", str(db), "--out", str(out)]) == 0
    assert pd.read_csv(out, sep="\t")["call"].tolist() == ["same"]  # LD came from the registered panel


def _study(world, causal_cols, seed):
    from scipy import stats

    from plocust.io import standardize_sumstats

    geno = Genotypes.open(world.bfile)
    v = geno.variants
    g = geno.read(v["idx"].to_numpy())
    rng = np.random.default_rng(seed)
    y = sum(g[:, c] / 2 for c in causal_cols) + rng.normal(0, 1, len(g))
    rows = []
    for j in range(g.shape[1]):
        r = stats.linregress(g[:, j], y) if g[:, j].std() > 0 else None
        rows.append(("1", int(v["pos"].iat[j]), v["id"].iat[j], v["a2"].iat[j], v["a1"].iat[j], v["a1"].iat[j],
                     r.slope if r else 0.0, r.stderr if r else np.nan, max(r.pvalue, 1e-300) if r else 1.0, len(g)))
    ss = standardize_sumstats(pd.DataFrame(rows, columns=["chrom", "pos", "id", "ref", "alt", "effect_allele",
                                                          "beta", "se", "p", "n"]))
    study = StudyInfo(species="Oryza sativa", trait="t", study=f"s{seed}", build="SRC", n_samples=len(g))
    return identify_loci(ss, study, genome=Genome(world.genome, "SRC"), geno=geno, cfg=IdentifyConfig(p_threshold=1e-6))


def test_coloc_same_and_distinct(world):
    geno = Genotypes.open(world.bfile)
    g = geno.read(geno.variants["idx"].to_numpy())
    c1 = int(np.flatnonzero(geno.variants["pos"] == CAUSAL_POS)[0])
    nxt = range((c1 // BLOCK_SNPS + 1) * BLOCK_SNPS, (c1 // BLOCK_SNPS + 2) * BLOCK_SNPS)
    c2 = max(nxt, key=lambda j: g[:, j].std())
    a, a2, b = _study(world, [c1], 1), _study(world, [c1], 2), _study(world, [c2], 3)
    same = compare(a, a2, "SRC", geno).sort_values("distance_bp").iloc[0]
    assert same["coloc_call"] == "same" and same["PP.H4"] > 0.8
    diff = compare(a, b, "SRC", geno).sort_values("distance_bp").iloc[0]
    assert diff["coloc_call"] == "distinct" and diff["PP.H3"] > 0.8


def test_coloc_across_builds(world, tmp_path):
    """Both studies moved to the target build by anchors; LD from a panel in target coordinates."""
    import shutil

    geno = Genotypes.open(world.bfile)
    c1 = int(np.flatnonzero(geno.variants["pos"] == CAUSAL_POS)[0])
    a, a2 = _study(world, [c1], 1), _study(world, [c1], 2)
    aligner = Aligner(world.target, build="TGT")
    anchor_passports(a, aligner)
    anchor_passports(a2, aligner)
    # the same panel, re-coordinated to the target build (5 kb insertion at 50 kb)
    bim = pd.read_csv(f"{world.bfile}.bim", sep="\t", header=None)
    bim[3] = np.where(bim[3] >= 50_000, bim[3] + INSERT_LEN, bim[3])
    bim.to_csv(tmp_path / "tgt.bim", sep="\t", header=False, index=False)
    for ext in ("bed", "fam"):
        shutil.copy(f"{world.bfile}.{ext}", tmp_path / f"tgt.{ext}")
    row = compare(a, a2, "TGT", Genotypes.open(tmp_path / "tgt")).sort_values("distance_bp").iloc[0]
    assert row["coloc_call"] == "same" and row["coloc_snps"] > 50


def test_imprint_territory(world):
    """Two loci: each imprint's territory stops halfway to the other lead."""
    geno = Genotypes.open(world.bfile)
    g = geno.read(geno.variants["idx"].to_numpy())
    c1 = int(np.flatnonzero(geno.variants["pos"] == CAUSAL_POS)[0])
    nxt = range((c1 // BLOCK_SNPS + 1) * BLOCK_SNPS, (c1 // BLOCK_SNPS + 2) * BLOCK_SNPS)
    c2 = max(nxt, key=lambda j: g[:, j].std())
    ps = sorted(_study(world, [c1, c2], 5), key=lambda p: p.signal.lead.pos)
    assert len(ps) >= 2
    a, b = ps[0], ps[1]
    mid = (a.signal.lead.pos + b.signal.lead.pos) // 2
    assert a.imprint.territory_end == mid and b.imprint.territory_start == mid + 1


def test_missing_lead_uses_best_associated_proxy(world, loci, tmp_path):
    from .conftest import write_bed

    geno = Genotypes.open(world.bfile)
    v = geno.variants
    lead = loci[0].signal.lead.pos
    keep = v[v["pos"] != lead]
    g = geno.read(keep["idx"].to_numpy())
    bim = keep.assign(cm=0)[["chrom", "id", "cm", "pos", "a1", "a2"]]
    write_bed(tmp_path / "nolead", g, bim, geno.samples["iid"].tolist())
    panel = Genotypes.open(tmp_path / "nolead")
    assert panel.nearest("1", lead, 0) is None
    row = compare([loci[0]], [loci[0].model_copy(deep=True)], "SRC", panel).iloc[0]
    assert row["call"] == "same" and row["r2"] > 0.9


def test_passports_from_hits(world):
    """A published hit table (no summary statistics, no genotypes) becomes anchored, matchable passports."""
    from plocust.hits import passports_from_hits
    from plocust.passport import LocusKind

    genome = Genome(world.genome, "SRC")
    genes = read_genes(world.gff)
    hits = pd.DataFrame({"chrom": ["chr01", "chr01"], "pos": [CAUSAL_POS, 200_500], "p": [1e-12, 1e-8],
                         "id": ["a", "b"], "gene": ["G1", "G3"]})
    ps = passports_from_hits(hits, genome, "Oryza sativa", "toy eQTL", kind=LocusKind.eqtl, genes=genes)
    assert [p.kind for p in ps] == [LocusKind.eqtl] * 2
    assert ps[0].trait.name == "expression of G1" and ps[0].genes[0].distance_bp == 0
    assert ps[0].lead_anchor.sequence[ps[0].lead_anchor.variant_offset] == genome.fetch("1", CAUSAL_POS, CAUSAL_POS)
    pl = anchor_passports(ps, Aligner(world.target, build="TGT"))[0]
    assert pl.lead_pos == CAUSAL_POS + INSERT_LEN


def test_merge_fragmented_clumps(world):
    """With a tiny clumping window one LD block splits into several clumps; merging restores one locus."""
    ss = read_sumstats(world.sumstats_path)
    geno = Genotypes.open(world.bfile)
    split = identify_loci(ss, STUDY, genome=Genome(world.genome, "SRC"), geno=geno,
                          cfg=IdentifyConfig(p_threshold=1e-8, window_kb=2, merge_r2=None))
    merged = identify_loci(ss, STUDY, genome=Genome(world.genome, "SRC"), geno=geno,
                           cfg=IdentifyConfig(p_threshold=1e-8, window_kb=2))
    assert len(split) > 1 and len(merged) == 1
    assert merged[0].signal.n_significant == sum(p.signal.n_significant for p in split)


def test_to_panel_exact_and_near(world):
    from plocust.coloc import to_panel

    geno = Genotypes.open(world.bfile)
    v = geno.variants.iloc[:3]
    df = pd.DataFrame({"chrom": "1", "pos": [v["pos"].iat[0], v["pos"].iat[1] + 7, 999_999],
                       "ref": [v["a2"].iat[0], v["a1"].iat[1], "A"], "alt": [v["a1"].iat[0], v["a2"].iat[1], "C"],
                       "z": [2.0, 3.0, 1.0]})
    m = to_panel(df, geno, tol=20).sort_values("pos")
    assert m["idx"].tolist() == [int(v["idx"].iat[0]), int(v["idx"].iat[1])]
    assert m["z_a1"].tolist() == [2.0, -3.0]  # second SNP reported for the other allele: sign flipped


def test_one_locus_many_traits_keeps_every_record(world, loci, tmp_path):
    """Same lead SNP for two traits: one passport_id (one locus), two database records (pleiotropy)."""
    a = loci[0]
    b = a.model_copy(deep=True)
    b.trait.name = "grain weight"
    assert a.passport_id == b.passport_id
    db = tmp_path / "db.sqlite"
    assert create(db, [a, b], name="t", version="0") == 2
    recs = LocusDB(db).records(a.passport_id)
    assert sorted(r.trait.name for r in recs) == ["grain weight", "plant height"]


def test_ld_call_without_coloc(world, loci):
    """The LD call must not depend on whether colocalization is switched on."""
    geno = Genotypes.open(world.bfile)
    other = loci[0].model_copy(deep=True)
    other.placements[0].lead_pos += SNP_STEP * BLOCK_SNPS
    other.credible_set = None
    for coloc in (True, False):
        rows = compare([loci[0]], [loci[0].model_copy(deep=True), other], "SRC", geno, MatchConfig(coloc=coloc))
        assert rows.sort_values("distance_bp")["call"].tolist() == ["same", "distinct_nearby"]


def test_finemap_passports_and_cli(world, loci, tmp_path):
    from plocust.finemap import finemap_passports
    from plocust.passport import write_passports

    ss = read_sumstats(world.sumstats_path)
    ps = [p.model_copy(deep=True) for p in loci]
    rows = finemap_passports(ps, ss, Genotypes.open(world.bfile), "SRC", n=300)
    fm = ps[0].fine_mapping
    assert rows and fm.signals and fm.lead_signal == 0
    assert any(v.pos == CAUSAL_POS for v in fm.signals[0].variants)  # the causal SNP is in the credible set
    f = tmp_path / "l.jsonl"
    write_passports(f, loci)
    out = tmp_path / "fm.jsonl"
    assert main(["finemap", str(f), "--sumstats", str(world.sumstats_path), "--bfile", str(world.bfile),
                 "--n", "300", "--out", str(out)]) == 0
    assert read_passports(out)[0].fine_mapping.signals


def test_loci_regions_merge_overlaps(loci):
    from plocust.finemap import loci_regions

    a = loci[0]
    b = a.model_copy(deep=True)
    b.placements[0].lead_pos += 100_000
    c = a.model_copy(deep=True)
    c.placements[0].lead_pos += 5_000_000
    regs = loci_regions([a, b, c], "SRC")
    assert [len(r[3]) for r in regs] == [2, 1]
