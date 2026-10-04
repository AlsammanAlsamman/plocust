"""Atlas step 2: the 44,354 leaf eQTL hits of Liu et al. 2022 (287 accessions) as passports (kind eqtl).

Published hits only: lead SNP + p-value + target gene. Positions are on IRGSP-1.0 (MSU7 coordinates).
Target genes are MSU loci, placed with the MSU7 locus table.
"""

from pathlib import Path

import pandas as pd

from plocust.anchor import Aligner, check_uniqueness
from plocust.hits import passports_from_hits
from plocust.io import Genome, normalize_chrom
from plocust.passport import LocusKind, write_passports

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = DATA / "atlas" / "passports"
RES = ROOT / "atlas" / "results"


def msu_genes() -> pd.DataFrame:
    t = pd.read_csv(DATA / "ref/msu7/locus_brief_info.7.0", sep="\t")
    g = t.groupby("locus").agg(chrom=("chr", "first"), start=("start", "min"), end=("stop", "max"),
                               annotation=("annotation", "first")).reset_index().rename(columns={"locus": "id"})
    g["chrom"] = g["chrom"].map(normalize_chrom)
    g["name"] = None
    return g


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    e = pd.read_csv(DATA / "expr/liu2022/liu2022_eqtl.tsv", sep="\t")
    hits = pd.DataFrame({"chrom": e["chrom"], "pos": e["pos"], "p": e["p"].astype(float).clip(lower=1e-300),
                         "id": e["snp"], "gene": e["gene"], "type": e["type"]})
    genome = Genome(DATA / "ref/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa", "IRGSP-1.0")
    ps = passports_from_hits(hits, genome, "Oryza sativa", "Liu et al. 2022 leaf eQTL", kind=LocusKind.eqtl,
                             genes=msu_genes(), n_samples=287)
    check_uniqueness(ps, Aligner(genome.path, build="IRGSP-1.0"))
    write_passports(OUT / "eqtl_liu2022.jsonl", ps)

    local = hits["type"].eq("L").to_numpy()
    unique = pd.Series([p.lead_anchor.unique for p in ps])
    dist = pd.Series([p.genes[0].distance_bp if p.genes else None for p in ps])
    s = pd.DataFrame([{
        "eqtl_hits": len(hits), "passports": len(ps), "unique_lead_anchor": round(unique.mean(), 3),
        "target_gene_placed": round(dist.notna().mean(), 3),
        "local_within_100kb_of_gene": round((dist[local] <= 100_000).mean(), 3),
    }])
    RES.mkdir(parents=True, exist_ok=True)
    s.to_csv(RES / "eqtl_passports.tsv", sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
