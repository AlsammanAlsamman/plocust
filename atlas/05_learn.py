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


def msu_to_rap(rap_genes: pd.DataFrame) -> dict:
    """MSU locus -> RAP gene with the largest overlap on IRGSP-1.0 (both annotations use that assembly)."""
    t = pd.read_csv(DATA / "ref/msu7/locus_brief_info.7.0", sep="\t")
    msu = t.groupby("locus").agg(chrom=("chr", "first"), start=("start", "min"), end=("stop", "max")).reset_index()
    msu["chrom"] = msu["chrom"].map(normalize_chrom)
    out = {}
    for chrom, m in msu.groupby("chrom"):
        r = rap_genes[rap_genes["chrom"] == chrom].sort_values("start")
        rs, re_, rid = r["start"].to_numpy(), r["end"].to_numpy(), r["id"].to_numpy()
        for locus, a, b in zip(m["locus"], m["start"], m["end"]):
            i = np.searchsorted(re_, a)  # candidates whose end >= a (ends are nearly sorted with starts)
            best, best_ov = None, 0
            for j in range(max(0, i - 5), min(len(rs), i + 10)):
                ov = min(b, re_[j]) - max(a, rs[j])
                if ov > best_ov:
                    best, best_ov = rid[j], ov
            if best is not None:
                out[locus] = best
    return out


def nearest_gene(genes: pd.DataFrame, chrom: str, pos: int):
    g = genes[(genes["chrom"] == chrom) & (genes["biotype"] == "protein_coding")]
    d = np.maximum(0, np.maximum(g["start"].to_numpy() - pos, pos - g["end"].to_numpy()))
    return g["id"].iat[int(np.argmin(d))] if len(g) else None


def d_eqtl(gwas, eqtl, dense, known, rap_genes):
    """Same-signal local eQTLs per GWAS locus, and the benchmark: eQTL gene vs nearest gene, against cloned genes."""
    m2r = msu_to_rap(rap_genes)
    idx = defaultdict(list)
    for e in eqtl:
        if e.genes and e.genes[0].distance_bp <= 100_000:  # local (cis) eQTLs only
            idx[e.placement(B).chrom].append(e)
    kn = {(r.panel, r.trait, r.passport_id): r.nearest_known_id for r in known.itertuples()} if len(known) else {}
    cfg = MatchConfig(max_distance_kb=250, proxy_kb=20, coloc=False)
    rows = []
    for panel, trait, p in gwas:
        pl = p.placement(B)
        if pl is None:
            continue
        near = [e for e in idx[pl.chrom] if abs(e.placement(B).lead_pos - pl.lead_pos) <= 250_000]
        gene_of = {e.passport_id: m2r.get(e.genes[0].id, e.genes[0].id) for e in near}
        same = []
        if near:
            pairs = compare([p], near, B, dense, cfg)
            hit = pairs[pairs["call"].isin(["same", "same_opposite_effect"])].sort_values("r2", ascending=False)
            same = list(dict.fromkeys(hit["target"].map(gene_of)))
        rows.append({"panel": panel, "trait": trait, "passport_id": p.passport_id, "eqtls_within_250kb": len(near),
                     "same_signal_eqtl_genes": ";".join(same[:5]), "n_same_genes": len(same),
                     "nearest_gene": nearest_gene(rap_genes, pl.chrom, pl.lead_pos),
                     "known_gene": kn.get((panel, trait, p.passport_id))})
    return pd.DataFrame(rows)


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

    from plocust.io import read_genes

    rap_genes = read_genes(DATA / "ref/Oryza_sativa.IRGSP-1.0.63.gff3.gz")
    d = d_eqtl(gwas, eqtl, dense, known, rap_genes)
    d.to_csv(RES / "learn_D_gwas_eqtl.tsv", sep="\t", index=False)
    with_eqtl = d[d["n_same_genes"] > 0]
    bench = with_eqtl[with_eqtl["known_gene"].notna()]
    top_eqtl = bench.apply(lambda r: r["same_signal_eqtl_genes"].split(";")[0] == r["known_gene"], axis=1)
    any_eqtl = bench.apply(lambda r: r["known_gene"] in r["same_signal_eqtl_genes"].split(";"), axis=1)
    nearest = bench["nearest_gene"] == bench["known_gene"]
    summary = pd.DataFrame([{
        "gwas_loci": len(d), "with_same_signal_eqtl": len(with_eqtl),
        "benchmark_loci (cloned gene + same-signal eQTL)": len(bench),
        "top eQTL gene = cloned gene": int(top_eqtl.sum()), "any eQTL gene = cloned gene": int(any_eqtl.sum()),
        "nearest gene = cloned gene": int(nearest.sum()),
        "loci without cloned gene but with eQTL gene (new candidates)": int(with_eqtl["known_gene"].isna().sum()),
    }]).T
    summary.to_csv(RES / "learn_D_benchmark.tsv", sep="\t", header=False)
    print("\nD. GWAS x eQTL\n", summary.to_string(header=False), flush=True)


if __name__ == "__main__":
    main()
