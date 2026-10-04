"""Atlas step 4: the rice passport database (gene + GWAS + eQTL passports, chip markers, LD panel)."""

from pathlib import Path

import pandas as pd

from plocust.db import LocusDB, add_markers, attach_ld_panel, create
from plocust.passport import read_passports

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PASS = DATA / "atlas" / "passports"
OUT = DATA / "db" / "plocust-db-rice-atlas.sqlite"
VERSION = "0.1.0-atlas"


def main():
    genes = LocusDB(DATA / "db/plocust-db-rice.sqlite")
    passports = [p for p in genes.by_trait("")]
    sources = {"gene (funRiceGenes)": len(passports)}
    for f in sorted(PASS.glob("gwas_*.jsonl")):
        ps = read_passports(f)
        passports += ps
        sources["gwas"] = sources.get("gwas", 0) + len(ps)
    eq = read_passports(PASS / "eqtl_liu2022.jsonl")
    passports += eq
    sources["eqtl (Liu et al. 2022)"] = len(eq)
    n = create(OUT, passports, name="plocust-db-rice", version=VERSION,
               extra_meta={"species": "Oryza sativa", "reference_build": "IRGSP-1.0",
                           "sources": "; ".join(f"{k}: {v}" for k, v in sources.items())})
    f7 = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU7.SNP_flanking_seq.txt", sep="\t", dtype=str)
    f7["sequence"] = f7["msu7_seq1"] + "[" + f7["alleles"] + "]" + f7["msu7_seq2"]
    m = add_markers(OUT, f7.rename(columns={"snp_id": "marker_id"}), "RDP1 44K array")
    attach_ld_panel(OUT, DATA / "3k/core3k", "IRGSP-1.0", "3K RG core SNPs (365K, 3,000 accessions)")
    info = LocusDB(OUT).info()
    print(f"{n} passports ({sources}); {m} markers; {OUT.stat().st_size / 1e6:.0f} MB")
    for k, v in info.items():
        print(f"  {k:24} {v}")


if __name__ == "__main__":
    main()
