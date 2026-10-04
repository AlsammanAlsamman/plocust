"""Validation B: do gene passports built on IRGSP-1.0 land on the right gene in other assemblies?

Truth: Ensembl Compara one-to-one orthologs (4 cultivar assemblies) and MSU6 gene models (old build).
Baseline: reuse the IRGSP-1.0 coordinates unchanged on the other assembly ("same coordinates").
A placement is correct when the anchored gene interval overlaps the true gene (+-2 kb).
"""

import gzip
import re
from pathlib import Path

import pandas as pd

from plocust.db import LocusDB
from plocust.io import normalize_chrom, read_genes

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "validation" / "results" / "b_anchoring"
OUT.mkdir(parents=True, exist_ok=True)
PAD = 2_000

ASSEMBLIES = {
    "oryza_sativa_mh63": ("MH63RS2", "asm/Oryza_sativa_mh63.MH63RS2.63.gff3.gz"),
    "oryza_sativa_zs97": ("ZS97RS3", "asm/Oryza_sativa_zs97.ZS97RS3.63.gff3.gz"),
    "oryza_sativa_n22": ("N22RS2", "asm/Oryza_sativa_n22.OsN22RS2.63.gff3.gz"),
    "oryza_sativa_azucena": ("AzucenaRS1", "asm/Oryza_sativa_azucena.AzucenaRS1.63.gff3.gz"),
}


def msu6_genes() -> pd.DataFrame:
    rows = []
    with gzip.open(DATA / "msu6/all.gff3.gz", "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene" and (m := re.search(r"Alias=(LOC_Os\w+)", f[8])):
                rows.append((normalize_chrom(f[0]), int(f[3]), int(f[4]), m.group(1)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "id"])


def overlaps(chrom, start, end, g) -> bool:
    return chrom == g["chrom"] and start <= g["end"] + PAD and end >= g["start"] - PAD


def main():
    db = LocusDB(DATA / "db/plocust-db-rice.sqlite")
    by_gene = {p.genes[0].id: p for p in db.by_trait("") if p.genes}
    orth = pd.read_csv(DATA / "known/orthologs.tsv", sep="\t")
    orth = orth[orth["type"] == "ortholog_one2one"]
    rows = []
    for species, (build, gff) in ASSEMBLIES.items():
        genes = read_genes(DATA / gff).set_index("id")
        for _, o in orth[orth["species"] == species].iterrows():
            p = by_gene.get(o["rap_id"])
            if p is None or o["target_gene"] not in genes.index:
                continue
            truth = genes.loc[o["target_gene"]]
            src, pl = p.placement("IRGSP-1.0"), p.placement(build)
            rows.append({"build": build, "gene": o["rap_id"], "placed": pl is not None,
                         "correct": pl is not None and overlaps(pl.chrom, pl.start, pl.end, truth),
                         "same_coords_correct": overlaps(src.chrom, src.start, src.end, truth),
                         "flags": ",".join(f.value for f in pl.flags) if pl else "unplaced"})
    # old build: MSU6 gene models, matched through the MSU locus IDs of funRiceGenes
    info = pd.read_csv(DATA / "known/geneInfo.table.txt", sep="\t", quoting=3, encoding_errors="replace", dtype=str)
    rap2msu = dict(zip(info["RAPdb"], info["MSU"]))
    m6 = msu6_genes().drop_duplicates("id").set_index("id")
    for rap in orth["rap_id"].unique():
        p, msu = by_gene.get(rap), rap2msu.get(rap)
        if p is None or not isinstance(msu, str) or msu not in m6.index:
            continue
        truth = m6.loc[msu]
        src, pl = p.placement("IRGSP-1.0"), p.placement("MSU6")
        rows.append({"build": "MSU6", "gene": rap, "placed": pl is not None,
                     "correct": pl is not None and overlaps(pl.chrom, pl.start, pl.end, truth),
                     "same_coords_correct": overlaps(src.chrom, src.start, src.end, truth),
                     "flags": ",".join(f.value for f in pl.flags) if pl else "unplaced"})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "gene_anchoring.tsv", sep="\t", index=False)
    s = df.groupby("build").agg(genes=("gene", "size"), placed=("placed", "mean"), correct=("correct", "mean"),
                                same_coords_correct=("same_coords_correct", "mean"))
    s["correct_of_placed"] = df[df["placed"]].groupby("build")["correct"].mean()
    s = s.round(3)
    s.to_csv(OUT / "summary.tsv", sep="\t")
    print(s.to_string())
    wrong = df[df["placed"] & ~df["correct"]]
    print(f"\nwrong placements: {len(wrong)}; flags among them: {wrong['flags'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
