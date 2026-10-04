"""Atlas step 3: passports for every GWAS trait (3K on IRGSP-1.0; RDP1 on MSU6, then anchored to IRGSP-1.0).

Thresholds: 3K Bonferroni with >= 2 significant SNPs per locus; RDP1 p < 1e-5 (44K array, ~400 lines).

RDP1 uses the published EMMA p-values: REGENIE over-corrects in this small, structured panel (lambda 0.37-0.72).
Effect signs (direction only) come from the REGENIE re-analysis.
"""

import re
import sys
import time
from pathlib import Path

import pandas as pd

from plocust.anchor import Aligner, anchor_passports, check_uniqueness
from plocust.identify import IdentifyConfig, StudyInfo, identify_loci
from plocust.io import Genome, Genotypes, read_genes, read_sumstats
from plocust.passport import write_passports

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
GWAS = DATA / "atlas" / "gwas"
OUT = DATA / "atlas" / "passports"
RES = ROOT / "atlas" / "results"

PANELS = {
    "3K": dict(build="IRGSP-1.0", genome="ref/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa", geno="3k/core3k",
               study="3K RG (IRGCIS phenotypes, REGENIE)", cfg=IdentifyConfig(min_significant=2)),
    "RDP1": dict(build="MSU6", genome="msu6/MSU6.fa", geno="rdp1/rdp1",
                 study="RDP1 (Zhao et al. 2011 phenotypes, REGENIE)", cfg=IdentifyConfig(p_threshold=1e-5)),
}


_PUBLISHED = None


def rdp1_published(regenie: pd.DataFrame, column: str) -> pd.DataFrame:
    """Published EMMA p-values for one RDP1 trait, with the effect sign of the REGENIE re-analysis."""
    global _PUBLISHED
    if _PUBLISHED is None:
        _PUBLISHED = pd.read_csv(DATA / "rdp1/MixedModel_Pval_all.txt", sep="\t")
        _PUBLISHED.columns = [c.strip() for c in _PUBLISHED.columns]
    key = re.sub(r"[^a-z0-9]", "", column.lower().replace("arberdeen", "aberdeen"))
    match = [c for c in _PUBLISHED.columns if re.sub(r"[^a-z0-9]", "", c.lower().replace("arberdeen", "aberdeen")) == key]
    if not match:
        print(f"  no published p-values for {column!r}; using REGENIE", flush=True)
        return regenie
    pub = _PUBLISHED[["SNPID", match[0]]].rename(columns={"SNPID": "id", match[0]: "p_pub"})
    m = regenie.merge(pub, on="id", how="inner")
    m["p"] = pd.to_numeric(m["p_pub"], errors="coerce")
    m["se"] = float("nan")  # REGENIE beta kept for its sign only
    return m.drop(columns="p_pub").dropna(subset=["p"]).query("p > 0").reset_index(drop=True)


def main(panels=("3K", "RDP1")):
    OUT.mkdir(parents=True, exist_ok=True)
    traits = pd.read_csv(RES / "traits.tsv", sep="\t")
    irgsp = Genome(DATA / PANELS["3K"]["genome"], "IRGSP-1.0")
    genes = read_genes(DATA / "ref/Oryza_sativa.IRGSP-1.0.63.gff3.gz")
    a_irgsp = Aligner(irgsp.path, build="IRGSP-1.0")
    rows = []
    for panel in panels:
        cfg = PANELS[panel]
        genome = Genome(DATA / cfg["genome"], cfg["build"])
        geno = Genotypes.open(DATA / cfg["geno"])
        a_src = a_irgsp if cfg["build"] == "IRGSP-1.0" else Aligner(genome.path, build=cfg["build"])
        for t in traits[traits["panel"] == panel].itertuples():
            f = GWAS / panel / f"gwas_{t.trait}.regenie"
            out = OUT / f"gwas_{panel}_{t.trait}.jsonl"
            if not f.exists():
                continue
            t0 = time.time()
            ss = read_sumstats(f, sep=" ")
            if panel == "RDP1":
                ss = rdp1_published(ss, t.source_column)
            study = StudyInfo("Oryza sativa", t.trait.replace("_", " "), cfg["study"], cfg["build"], n_samples=int(t.n))
            ps = identify_loci(ss, study, genome=genome, genes=genes if cfg["build"] == "IRGSP-1.0" else None,
                               geno=geno, cfg=cfg["cfg"])
            check_uniqueness(ps, a_src)
            if cfg["build"] != "IRGSP-1.0":
                anchor_passports(ps, a_irgsp)
            write_passports(out, ps)
            placed = sum(p.placement("IRGSP-1.0") is not None for p in ps)
            rows.append({"panel": panel, "trait": t.trait, "loci": len(ps), "placed_on_IRGSP": placed,
                         "unique_lead_anchor": sum(bool(p.lead_anchor.unique) for p in ps),
                         "median_block_kb": round(pd.Series([p.ld_block.length_kb for p in ps if p.ld_block]).median(), 1)
                         if ps else None,
                         "seconds": round(time.time() - t0, 1)})
            print(rows[-1], flush=True)
    s = pd.DataFrame(rows)
    s.to_csv(RES / f"gwas_passports_{'_'.join(panels)}.tsv", sep="\t", index=False)
    print(s.groupby("panel")[["loci", "placed_on_IRGSP", "unique_lead_anchor", "seconds"]].sum().to_string())


if __name__ == "__main__":
    main(tuple(sys.argv[1:]) or ("3K", "RDP1"))
