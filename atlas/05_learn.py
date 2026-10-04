"""Atlas step 5: run the tool on everything and learn.

  A  loci per trait, and fragmentation (separate loci of one trait whose leads are in LD)
  B  known genes: cloned genes whose literature matches the trait, inside each locus (+-50 kb)
  C  cross-study: 3K vs RDP1 for paired traits (same-signal test, dense 3K LD)
  D  GWAS x eQTL: same-signal eQTLs at each GWAS locus -> candidate genes;
     benchmark at loci with a trait-matched cloned gene: eQTL gene vs nearest gene
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from plocust.db import LocusDB
from plocust.io import Genotypes, normalize_chrom
from plocust.ld import r2_matrix
from plocust.match import MatchConfig, compare
from plocust.passport import read_passports

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PASS = DATA / "atlas" / "passports"
RES = ROOT / "atlas" / "results"
B = "IRGSP-1.0"

# trait -> regex over funRiceGenes keywords + paper titles
TRAIT_TERMS = {
    "culm_length|plant_height|seedling_height": r"plant height|dwarf|semi-dwarf|culm length|internode|gibberellin",
    "glutinous_endosperm|amylose_content|alkali_spreading_value": r"amylose|waxy|glutinous|eating|cooking|gelatini|starch",
    "seed_coat_colour|pericarp_color": r"pericarp|proanthocyanidin|red rice|seed colou?r",
    r"apiculus|hull_colour|leaf_sheath|auricle|collar|sterile_lemma_colour|internode_colour|ligule|node_colour":
        r"anthocyanin|pigment|colou?r|purple|flavonoid",
    "pubescence": r"trichome|pubescen|hairy|hair|glabrous",
    "awn": r"\bawn",
    "shattering|threshability": r"shatter|abscission",
    "culm_angle|leaf_angle|flag_leaf_angle|culm_habit": r"tiller angle|leaf angle|lamina|prostrate|erect|gravitrop|brassinosteroid",
    "panicle|branching|exsertion|florets|seed_number": r"panicle|inflorescence|branch|spikelet",
    "flowering|heading": r"heading|flowering|photoperiod",
    "seed_length|seed_width|seed_volume|seed_surface|brown_rice|length_width": r"grain size|grain length|grain width|seed size|grain shape|grain weight",
    "blast": r"blast|magnaporthe|pyricularia|disease resistance",
    "protein_content": r"protein content|grain protein|storage protein",
    "straighthead": r"straighthead|arsenic",
    "culm_number|panicle_number": r"tiller",
    "leaf_length|flag_leaf|leaf_width": r"leaf size|flag leaf|leaf width|leaf length|narrow leaf",
}
PAIRS = {  # 3K trait -> RDP1 trait
    "culm_length": "plant_height", "glutinous_endosperm": "amylose_content", "seed_coat_colour": "pericarp_color",
    "awn": "awn_presence", "panicle_length": "panicle_length", "lemma_pubescence": "leaf_pubescence",
    "blade_pubescence": "leaf_pubescence", "culm_angle": "culm_habit", "panicle_shattering": "panicle_fertility",
}


def trait_regex(trait: str):
    for key, rx in TRAIT_TERMS.items():
        if re.search(key, trait):
            return re.compile(rx, re.I)
    return None


def gene_text() -> dict[str, str]:
    k = pd.read_csv(DATA / "known/geneKeyword.table.txt", sep="\t", quoting=csv.QUOTE_NONE,
                    encoding_errors="replace", dtype=str)
    txt = defaultdict(str)
    for gid, kw, title in zip(k["RAPdb"], k["Keyword"], k["Title"]):
        if isinstance(gid, str):
            txt[gid] += f" {kw} {title}"
    return txt


def load(pattern):
    ps = []
    for f in sorted(PASS.glob(pattern)):
        panel, trait = f.stem.split("_", 2)[1:]
        for p in read_passports(f):
            ps.append((panel, trait, p))
    return ps


def a_loci(gwas, core):
    rows = []
    by = defaultdict(list)
    for panel, trait, p in gwas:
        by[(panel, trait)].append(p)
    for (panel, trait), ps in by.items():
        frag = pairs = 0
        placed = [p for p in ps if p.placement(B)]
        for i, a in enumerate(placed):
            for b in placed[i + 1:]:
                pa, pb = a.placement(B), b.placement(B)
                if pa.chrom != pb.chrom or abs(pa.lead_pos - pb.lead_pos) > 1_000_000:
                    continue
                va, vb = core.nearest(pa.chrom, pa.lead_pos, 2000), core.nearest(pb.chrom, pb.lead_pos, 2000)
                if va is None or vb is None:
                    continue
                g = core.read([int(va["idx"]), int(vb["idx"])])
                pairs += 1
                frag += r2_matrix(g[:, :1], g[:, 1:])[0, 0] >= 0.5
        rows.append({"panel": panel, "trait": trait, "loci": len(ps), "pairs_within_1Mb": pairs,
                     "pairs_in_LD_r2>=0.5": int(frag)})
    return pd.DataFrame(rows)


def b_known(gwas, db, text):
    rows = []
    for panel, trait, p in gwas:
        rx = trait_regex(trait)
        pl = p.placement(B)
        if rx is None or pl is None:
            continue
        hits = []
        for k in db.region(B, pl.chrom, pl.start - 50_000, pl.end + 50_000):
            if k.kind.value == "gene" and rx.search(text.get(k.genes[0].id, "")):
                kp = k.placement(B)
                hits.append((abs(kp.lead_pos - pl.lead_pos), k.genes[0].id, k.genes[0].name))
        hits.sort()
        rows.append({"panel": panel, "trait": trait, "passport_id": p.passport_id, "chrom": pl.chrom,
                     "lead": pl.lead_pos, "interval_kb": round((pl.end - pl.start) / 1000, 1), "p": p.signal.pvalue,
                     "known_genes": ";".join(f"{n}({g}, {d // 1000} kb)" for d, g, n in hits[:5]),
                     "n_known": len(hits), "nearest_known_id": hits[0][1] if hits else None,
                     "nearest_known_kb": round(hits[0][0] / 1000, 1) if hits else None})
    return pd.DataFrame(rows)


def c_cross(gwas, dense):
    by = defaultdict(list)
    for panel, trait, p in gwas:
        by[(panel, trait)].append(p)
    rows, allpairs = [], []
    for t3, tr in PAIRS.items():
        a, b = by.get(("RDP1", tr), []), by.get(("3K", t3), [])
        if not a or not b:
            rows.append({"3K": t3, "RDP1": tr, "rdp1_loci": len(a), "3k_loci": len(b)})
            continue
        pairs = compare(a, b, B, dense, MatchConfig(max_distance_kb=1000, proxy_kb=20))
        pairs.insert(0, "pair", f"{tr}~{t3}")
        allpairs.append(pairs)
        best = pairs.sort_values("score", ascending=False).groupby("query").head(1) if len(pairs) else pairs
        rows.append({"3K": t3, "RDP1": tr, "rdp1_loci": len(a), "3k_loci": len(b),
                     "rdp1_with_3k_within_1Mb": best["query"].nunique() if len(best) else 0,
                     **{f"best_{c}": int((best["call"] == c).sum()) for c in
                        ("same", "same_opposite_effect", "ambiguous", "distinct_nearby")},
                     "coloc_same": int((best["coloc_call"] == "same").sum()) if len(best) else 0})
    if allpairs:
        pd.concat(allpairs).to_csv(RES / "learn_C_pairs.tsv", sep="\t", index=False)
    return pd.DataFrame(rows)


def d_eqtl(gwas, eqtl, dense, known):
    """Same-signal eQTLs per GWAS locus (local eQTLs only), and the benchmark against cloned genes."""
    idx = defaultdict(list)
    for e in eqtl:
        pl = e.placement(B)
        if e.genes and e.genes[0].distance_bp <= 100_000:  # local (cis) eQTLs
            idx[pl.chrom].append(e)
    rows = []
    cfg = MatchConfig(max_distance_kb=250, proxy_kb=20, coloc=False)
    kn = known.set_index("passport_id") if len(known) else known
    for panel, trait, p in gwas:
        pl = p.placement(B)
        if pl is None:
            continue
        near = [e for e in idx[pl.chrom] if abs(e.placement(B).lead_pos - pl.lead_pos) <= 250_000]
        if not near:
            continue
        pairs = compare([p], near, B, dense, cfg)
        same = pairs[pairs["call"].isin(["same", "same_opposite_effect"])]
        gene_of = {e.passport_id: e.genes[0].id for e in near}
        same_genes = same.sort_values("r2", ascending=False)["target"].map(gene_of).drop_duplicates().tolist()
        k = kn.loc[p.passport_id] if len(kn) and p.passport_id in kn.index else None
        rows.append({"panel": panel, "trait": trait, "passport_id": p.passport_id, "eqtls_within_250kb": len(near),
                     "same_signal_eqtl_genes": ";".join(same_genes[:5]), "n_same_genes": len(same_genes),
                     "known_gene": None if k is None else k["nearest_known_id"]})
    return pd.DataFrame(rows)


def msu_to_rap() -> dict:
    g = pd.read_csv(DATA / "known/geneInfo.table.txt", sep="\t", quoting=csv.QUOTE_NONE, encoding_errors="replace", dtype=str)
    return {m: r for m, r in zip(g["MSU"], g["RAPdb"]) if isinstance(m, str) and isinstance(r, str)}


def main():
    RES.mkdir(parents=True, exist_ok=True)
    gwas = load("gwas_*.jsonl")
    eqtl = read_passports(PASS / "eqtl_liu2022.jsonl")
    core, dense = Genotypes.open(DATA / "3k/core3k"), Genotypes.open(DATA / "3k/filt/base_filtered_v0.7")
    db = LocusDB(DATA / "db/plocust-db-rice.sqlite")
    print(f"{len(gwas)} GWAS passports, {len(eqtl)} eQTL passports", flush=True)

    a = a_loci(gwas, core)
    a.to_csv(RES / "learn_A_loci.tsv", sep="\t", index=False)
    print("\nA. LOCI AND FRAGMENTATION\n", a[a["loci"] > 0].to_string(index=False), flush=True)

    known = b_known(gwas, db, gene_text())
    known.to_csv(RES / "learn_B_known_genes.tsv", sep="\t", index=False)
    kb = known.groupby(["panel", "trait"]).agg(loci=("passport_id", "size"), with_known_gene=("n_known", lambda x: int((x > 0).sum())))
    print("\nB. LOCI WITH A TRAIT-MATCHED CLONED GENE (+-50 kb)\n", kb.to_string(), flush=True)

    c = c_cross(gwas, dense)
    c.to_csv(RES / "learn_C_cross_study.tsv", sep="\t", index=False)
    print("\nC. 3K vs RDP1\n", c.to_string(index=False), flush=True)

    d = d_eqtl(gwas, eqtl, dense, known)
    rap = msu_to_rap()
    d["same_genes_rap"] = d["same_signal_eqtl_genes"].fillna("").map(lambda s: ";".join(rap.get(g, g) for g in s.split(";") if g))
    d.to_csv(RES / "learn_D_gwas_eqtl.tsv", sep="\t", index=False)
    bench = d[d["known_gene"].notna() & (d["n_same_genes"] > 0)]
    hit = bench.apply(lambda r: r["known_gene"] in r["same_genes_rap"].split(";"), axis=1) if len(bench) else pd.Series(dtype=bool)
    print(f"\nD. GWAS x eQTL: {len(d)} loci with local eQTLs nearby; {int((d['n_same_genes'] > 0).sum())} with a "
          f"same-signal eQTL gene; benchmark loci (known gene + same-signal eQTL): {len(bench)}, "
          f"eQTL gene = known gene: {int(hit.sum()) if len(bench) else 0}", flush=True)


if __name__ == "__main__":
    main()
