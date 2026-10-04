"""plocust: build-independent locus passports for plant GWAS."""

import argparse
import csv
import json
import sys
from pathlib import Path

import pandas as pd
from pydantic import ValidationError

from . import __version__
from .passport import LocusPassport, json_schema, read_passports, write_passports


def _kv(pairs):
    return dict(p.split("=", 1) for p in pairs or [])


def _summary(passports, build) -> pd.DataFrame:
    rows = []
    for p in passports:
        pl = p.placement(build) if build else p.placements[0]
        rows.append({
            "passport_id": p.passport_id, "kind": p.kind.value, "trait": p.trait.name,
            "lead": p.signal.lead.id if p.signal else (p.genes[0].name if p.genes else None),
            "p": p.signal.pvalue if p.signal else None,
            "build": pl.build if pl else build, "chrom": pl.chrom if pl else None,
            "start": pl.start if pl else None, "end": pl.end if pl else None, "lead_pos": pl.lead_pos if pl else None,
            "flags": ",".join(f.value for f in pl.flags) if pl else "unplaced",
            "credible_set": len(p.credible_set.variants) if p.credible_set else None,
            "block_kb": p.ld_block.length_kb if p.ld_block else None,
            "n_haplotypes": p.ld_block.n_haplotypes if p.ld_block else None,
            "nearest_gene": p.genes[0].id if p.genes else None,
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ commands


def cmd_schema(args):
    print(json.dumps(json_schema(), indent=2))
    return 0


def cmd_validate(args):
    status = 0
    for path in args.passports:
        try:
            ps = read_passports(path)
        except (ValidationError, ValueError, OSError) as e:
            print(f"FAIL  {path}\n{e}", file=sys.stderr)
            status = 1
        else:
            print(f"OK    {path}  {len(ps)} passport(s)" if len(ps) != 1 else f"OK    {path}  {ps[0].passport_id}")
    return status


def cmd_detect_build(args):
    from .build import detect_build
    from .io import Genome, read_sumstats

    ss = read_sumstats(args.sumstats, _kv(args.column), sep=args.sep)
    genomes = []
    for spec in args.genome:
        fasta, build = spec.rsplit(":", 1) if ":" in Path(spec).name else (spec, None)
        genomes.append(Genome(fasta, build=build))
    print(detect_build(ss, genomes, n=args.n).to_string(index=False))
    return 0


def cmd_identify(args):
    from .identify import IdentifyConfig, StudyInfo, identify_loci
    from .io import FlankTable, Genome, Genotypes, read_genes, read_sumstats

    ss = read_sumstats(args.sumstats, _kv(args.column), sep=args.sep)
    genome = Genome(args.genome, build=args.build) if args.genome else None
    build = args.build or (genome.build if genome else None)
    if build is None:
        sys.exit("--build is required when no --genome is given")
    flanks = FlankTable.read(args.flanks) if args.flanks else None
    if flanks is None and args.db:
        from .db import LocusDB

        flanks = LocusDB(args.db).marker_flanks()
    cfg = IdentifyConfig(p_threshold=args.p, window_kb=args.window_kb, flank=args.flank, block_r2=args.block_r2,
                         min_significant=args.min_significant)
    study = StudyInfo(species=args.species, trait=args.trait, study=args.study, build=build,
                      trait_ontology=args.ontology, panel=args.panel, n_samples=args.n, gwas_method=args.method,
                      doi=args.doi)
    passports = identify_loci(ss, study, genome=genome, flanks=flanks,
                              genes=read_genes(args.gff) if args.gff else None,
                              geno=Genotypes.open(args.bfile) if args.bfile else None, cfg=cfg)
    if args.check_unique and genome:
        from .anchor import Aligner, check_uniqueness

        check_uniqueness(passports, Aligner(args.genome, build=build))
    write_passports(args.out, passports)
    _summary(passports, build).to_csv(Path(args.out).with_suffix(".tsv"), sep="\t", index=False)
    print(f"{len(passports)} loci -> {args.out}")
    return 0


def cmd_anchor(args):
    from .anchor import Aligner, anchor_passports

    passports = read_passports(args.passports)
    aligner = Aligner(args.genome, build=args.build)
    placed = anchor_passports(passports, aligner, min_mapq=args.min_mapq)
    write_passports(args.out, passports)
    _summary(passports, aligner.build).to_csv(Path(args.out).with_suffix(".tsv"), sep="\t", index=False)
    print(f"placed {sum(p is not None for p in placed)}/{len(placed)} on {aligner.build} -> {args.out}")
    return 0


def _geno(args):
    """Genotypes for LD: --bfile, else the reference panel registered in --db for --build."""
    from .io import Genotypes

    if args.bfile:
        return Genotypes.open(args.bfile)
    if getattr(args, "db", None):
        from .db import LocusDB

        g = LocusDB(args.db).ld_panel(args.build)
        if g is not None:
            print(f"LD from the reference panel registered in {args.db}", file=sys.stderr)
        return g
    return None


def cmd_compare(args):
    from .match import MatchConfig, compare

    pairs = compare(read_passports(args.a), read_passports(args.b), args.build, _geno(args),
                    MatchConfig(max_distance_kb=args.max_kb))
    pairs.to_csv(args.out, sep="\t", index=False)
    print(pairs["call"].value_counts().to_string() if len(pairs) else "no pairs within range")
    return 0


def cmd_match(args):
    from .db import LocusDB
    from .match import MatchConfig, compare

    passports = read_passports(args.passports)
    db = LocusDB(args.db)
    known = db.near(passports, args.build, args.max_kb)
    pairs = compare(passports, known, args.build, _geno(args), MatchConfig(max_distance_kb=args.max_kb))
    sym = {k.passport_id: (k.genes[0].name if k.genes else None) for k in known}
    traits = {k.passport_id: "|".join([k.trait.name] + [t.name for t in k.traits_other]) for k in known}
    pairs["target_name"] = pairs["target"].map(sym)
    pairs["target_traits"] = pairs["target"].map(traits)
    pairs.to_csv(args.out, sep="\t", index=False)
    print(f"{len(pairs)} matches for {pairs['query'].nunique() if len(pairs) else 0}/{len(passports)} loci -> {args.out}")
    return 0


def cmd_db_build_genes(args):
    from .anchor import Aligner, anchor_passports
    from .db import create, gene_passports
    from .io import Genome, read_genes

    genome = Genome(args.genome, build=args.build)
    read = dict(sep="\t", quoting=csv.QUOTE_NONE, encoding_errors="replace", dtype=str)
    table = pd.read_csv(args.genes, **read)
    kw = pd.read_csv(args.keywords, **read) if args.keywords else None
    passports = gene_passports(read_genes(args.gff), table, genome, args.species, args.study, keywords=kw)
    for spec in args.assembly or []:
        fasta, build = spec.rsplit(":", 1) if ":" in spec else (spec, None)
        aligner = Aligner(fasta, build=build)
        placed = anchor_passports(passports, aligner)
        print(f"  {aligner.build}: placed {sum(p is not None for p in placed)}/{len(placed)}")
    n = create(args.out, passports, name=args.name, version=args.version,
               extra_meta={"species": args.species, "source": args.study, "reference_build": genome.build})
    print(f"{n} gene passports -> {args.out}")
    return 0


def cmd_db_add_markers(args):
    from .db import add_markers

    t = pd.read_csv(args.markers, sep=args.sep, dtype=str)
    n = add_markers(args.db, t, args.panel, id_col=args.id_col or t.columns[0],
                    seq_col=args.seq_col or next(c for c in t.columns if "seq" in c.lower()))
    print(f"{n} markers ({args.panel}) -> {args.db}")
    return 0


def cmd_db_attach_ld(args):
    from .db import attach_ld_panel

    attach_ld_panel(args.db, args.bfile, args.build, args.name)
    print(f"LD panel {args.name} ({args.build}) registered in {args.db}")
    return 0


def cmd_db_info(args):
    from .db import LocusDB

    for k, v in LocusDB(args.db).info().items():
        print(f"{k:16} {v}")
    return 0


def cmd_db_query(args):
    from .db import LocusDB

    db = LocusDB(args.db)
    if args.region:
        chrom, rng = args.region.split(":")
        start, end = (int(x.replace(",", "")) for x in rng.split("-"))
        found = db.region(args.build, chrom, start, end)
    else:
        found = db.by_trait(args.trait)
    print(_summary(found, args.build).to_csv(sep="\t", index=False), end="")
    return 0


def cmd_db_download(args):
    from .db import download

    print(download(args.name, args.dir))
    return 0


def cmd_card(args):
    from .db import LocusDB
    from .io import Genotypes, read_genes, read_sumstats
    from .plot import locus_card

    passports = read_passports(args.passports)
    p = next((x for x in passports if x.passport_id == args.id), None) if args.id else passports[args.index]
    if p is None:
        sys.exit(f"{args.id} not found in {args.passports}")
    known = None
    if args.db:
        build = args.build or p.placements[0].build
        known = [k for k in LocusDB(args.db).near([p], build, args.pad_kb) if k.kind.value == "gene"]
    locus_card(p, args.out, sumstats=read_sumstats(args.sumstats, _kv(args.column)) if args.sumstats else None,
               genes=read_genes(args.gff) if args.gff else None,
               geno=Genotypes.open(args.bfile) if args.bfile else None, known=known, build=args.build,
               p_threshold=args.p, pad_kb=args.pad_kb)
    print(args.out)
    return 0


# ------------------------------------------------------------------- parser


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="plocust", description=__doc__)
    ap.add_argument("--version", action="version", version=f"plocust {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    sub.add_parser("schema", help="print the passport JSON Schema").set_defaults(func=cmd_schema)

    p = sub.add_parser("validate", help="check passport files (.json / .jsonl)")
    p.add_argument("passports", nargs="+")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("detect-build", help="which genome build do the GWAS coordinates belong to?")
    p.add_argument("--sumstats", required=True)
    p.add_argument("--column", action="append", metavar="STD=FILECOL")
    p.add_argument("--sep", default=None)
    p.add_argument("--genome", action="append", required=True, metavar="FASTA[:BUILD]", help="candidate (repeat)")
    p.add_argument("--n", type=int, default=2000, help="SNPs to test")
    p.set_defaults(func=cmd_detect_build)

    p = sub.add_parser("identify", help="detect loci in GWAS results and write passports")
    p.add_argument("--sumstats", required=True)
    p.add_argument("--column", action="append", metavar="STD=FILECOL", help="map a column, e.g. p=MLM_P")
    p.add_argument("--sep", default=None, help="field separator (default: guess)")
    p.add_argument("--genome", help="reference FASTA (indexed) of the GWAS build")
    p.add_argument("--flanks", help="per-SNP flanking sequences (TSV: id, sequence with [A/G])")
    p.add_argument("--db", help="known-loci database: look up marker flanking sequences by SNP ID")
    p.add_argument("--build", help="build name (default: from the FASTA name)")
    p.add_argument("--gff", help="gene annotation (GFF3) of the same build")
    p.add_argument("--bfile", help="PLINK genotypes of the GWAS panel, for LD")
    p.add_argument("--p", type=float, default=None, help="significance threshold (default: Bonferroni)")
    p.add_argument("--window-kb", type=float, default=500)
    p.add_argument("--min-significant", type=int, default=1)
    p.add_argument("--flank", type=int, default=100)
    p.add_argument("--block-r2", type=float, default=0.5)
    p.add_argument("--species", required=True)
    p.add_argument("--trait", required=True)
    p.add_argument("--ontology")
    p.add_argument("--study", required=True)
    p.add_argument("--panel")
    p.add_argument("--n", type=int)
    p.add_argument("--method")
    p.add_argument("--doi")
    p.add_argument("--check-unique", action="store_true", help="map anchors back to the genome (needs --genome)")
    p.add_argument("--out", required=True, help="output .jsonl (a .tsv summary is written next to it)")
    p.set_defaults(func=cmd_identify)

    p = sub.add_parser("anchor", help="place passports on another assembly by sequence")
    p.add_argument("passports")
    p.add_argument("--genome", required=True, help="target assembly FASTA")
    p.add_argument("--build")
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_anchor)

    p = sub.add_parser("compare", help="match loci between two studies")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--build", required=True, help="build both studies are placed on")
    p.add_argument("--bfile", help="genotype panel on that build, for the same-signal test")
    p.add_argument("--db", help="use the reference LD panel registered in this database when --bfile is not given")
    p.add_argument("--max-kb", type=float, default=1000)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("match", help="match loci against the known-loci database")
    p.add_argument("passports")
    p.add_argument("--db", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--bfile")
    p.add_argument("--max-kb", type=float, default=500)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_match)

    db = sub.add_parser("db", help="known-loci database").add_subparsers(dest="db_command", required=True)
    p = db.add_parser("build-genes", help="build a database of known genes")
    p.add_argument("--genes", required=True, help="TSV with gene IDs and symbols (e.g. funRiceGenes geneInfo)")
    p.add_argument("--keywords", help="TSV of gene trait keywords")
    p.add_argument("--gff", required=True)
    p.add_argument("--genome", required=True)
    p.add_argument("--build")
    p.add_argument("--assembly", action="append", metavar="FASTA[:BUILD]", help="also place genes on this assembly")
    p.add_argument("--species", required=True)
    p.add_argument("--study", default="funRiceGenes")
    p.add_argument("--name", default="plocust-db")
    p.add_argument("--version", default="dev")
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_db_build_genes)
    p = db.add_parser("add-markers", help="add SNP-chip / KASP marker flanking sequences")
    p.add_argument("db")
    p.add_argument("--markers", required=True, help="table with marker ID and sequence written as ...[A/G]...")
    p.add_argument("--panel", required=True, help="marker panel name, e.g. 'RDP1 44K'")
    p.add_argument("--id-col")
    p.add_argument("--seq-col")
    p.add_argument("--sep", default="\t")
    p.set_defaults(func=cmd_db_add_markers)
    p = db.add_parser("attach-ld", help="register a reference genotype panel for LD")
    p.add_argument("db")
    p.add_argument("--bfile", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--name", required=True)
    p.set_defaults(func=cmd_db_attach_ld)
    p = db.add_parser("info")
    p.add_argument("db")
    p.set_defaults(func=cmd_db_info)
    p = db.add_parser("query")
    p.add_argument("db")
    p.add_argument("--build", required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--region", help="chr:start-end")
    g.add_argument("--trait", help="trait keyword")
    p.set_defaults(func=cmd_db_query)
    p = db.add_parser("download")
    p.add_argument("name")
    p.add_argument("--dir", default=".")
    p.set_defaults(func=cmd_db_download)

    p = sub.add_parser("card", help="draw the locus card for one passport")
    p.add_argument("passports")
    p.add_argument("--id")
    p.add_argument("--index", type=int, default=0)
    p.add_argument("--sumstats")
    p.add_argument("--column", action="append", metavar="STD=FILECOL")
    p.add_argument("--gff")
    p.add_argument("--bfile")
    p.add_argument("--db")
    p.add_argument("--build")
    p.add_argument("--p", type=float)
    p.add_argument("--pad-kb", type=float, default=150)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_card)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
