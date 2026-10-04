"""Atlas step 10: export the data for the rice atlas figure (docs/atlas_figure/atlas_data.js)."""

import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from plocust.db import LocusDB
from plocust.passport import read_passports

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RES = ROOT / "atlas" / "results"
OUT = ROOT / "docs" / "atlas_figure" / "atlas_data.js"
B = "IRGSP-1.0"
BIN = 2_000_000

ORGANS = {  # organ -> trait regex (trait names as in atlas/results/traits.tsv)
    "grain": r"glutinous|seed_coat|apiculus|hull|amylose|alkali|seed_length|seed_width|seed_volume|seed_surface|"
             r"brown_rice|length_width|pericarp|protein|sterile_lemma|lemma_pubescence",
    "panicle": r"panicle|branch|threshab|shatter|exsertion|awn|spikelet|florets|seed_number|flowering|ft_ratio",
    "leaf": r"leaf|blade|sheath|auricle|collar|ligule|senescence|blast",
    "culm": r"culm|plant_height|internode|seedling|straighthead|node",
}


def organ(trait: str) -> str:
    for o, rx in ORGANS.items():
        if re.search(rx, trait):
            return o
    return "culm"


def main():
    fai = pd.read_csv(DATA / "ref/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.fai", sep="\t", header=None)
    chroms = [{"name": str(c), "length": int(L)} for c, L in zip(fai[0], fai[1]) if str(c).isdigit()]
    chroms.sort(key=lambda c: int(c["name"]))

    db = LocusDB(DATA / "db/plocust-db-rice-atlas.sqlite")
    symbol = {p.genes[0].id: p.genes[0].name for p in db.by_trait("") if p.kind.value == "gene" and p.genes}
    known = pd.read_csv(RES / "learn_B_known_genes.tsv", sep="\t")
    known_of = {(r.panel, r.trait, r.passport_id): r.nearest_known_id for r in known.itertuples()
                if isinstance(r.nearest_known_id, str)}

    loci, by_pid = [], defaultdict(list)
    for f in sorted((DATA / "atlas/passports_fm").glob("gwas_*.jsonl")):
        panel, trait = f.stem.split("_", 2)[1:]
        for p in read_passports(f):
            pl = p.placement(B)
            if pl is None:
                continue
            fm = p.fine_mapping
            sig = fm.signals[fm.lead_signal] if fm and fm.lead_signal is not None else None
            gid = known_of.get((panel, trait, p.passport_id))
            rec = {"chrom": pl.chrom, "pos": pl.lead_pos, "trait": trait.replace("_", " "), "panel": panel,
                   "organ": organ(trait), "mlog10p": round(min(-np.log10(p.signal.pvalue), 60), 2),
                   "cs_size": len(sig.variants) if sig else None,
                   "top_pip": round(sig.variants[0].pip, 3) if sig else None,
                   "gene": (symbol.get(gid) or gid) if gid else None}
            loci.append(rec)
            by_pid[p.passport_id].append(len(loci) - 1)

    # pleiotropy: one lead SNP (same passport_id) for several traits; cross-study same-signal pairs
    links = []
    for idx in by_pid.values():
        traits = {loci[i]["trait"] for i in idx}
        if len(traits) > 1:
            links += [{"a": idx[0], "b": j, "kind": "same SNP"} for j in idx[1:] if loci[j]["trait"] != loci[idx[0]]["trait"]]
    pairs = pd.read_csv(RES / "learn_C_pairs.tsv", sep="\t")
    key = {(l["chrom"], l["pos"]): i for i, l in enumerate(loci)}
    for r in pairs[pairs["call"].isin(["same", "same_opposite_effect"])].itertuples():
        a, b = key.get((str(r.chrom), int(r.query_lead))), key.get((str(r.chrom), int(r.target_lead)))
        if a is not None and b is not None and a != b:
            links.append({"a": a, "b": b, "kind": "same signal (RDP1 ~ 3K)"})

    # expression: local leaf eQTLs per 2 Mb
    e = pd.read_csv(DATA / "expr/liu2022/liu2022_eqtl.tsv", sep="\t")
    e = e[e["type"] == "L"]
    eqtl = []
    for c in chroms:
        sub = e[e["chrom"].astype(str) == c["name"]]
        counts = np.bincount((sub["pos"] // BIN).astype(int), minlength=c["length"] // BIN + 1)
        eqtl.append({"chrom": c["name"], "counts": counts.tolist()})

    # landmark genes: classic cloned genes, with the organ of their trait, drawn at their true positions
    landmarks_wanted = {"sd1": "culm", "Wx": "grain", "Rc": "grain", "GS3": "grain", "OsC1": "grain", "HL6": "leaf",
                        "ALK": "grain", "Pita2": "leaf", "GSE5": "grain", "Hd1": "panicle", "Ghd7": "panicle",
                        "GW5": "grain", "An-1": "panicle", "qSH1": "panicle", "PROG1": "culm"}
    landmarks = []
    for p in db.by_trait(""):
        if p.kind.value != "gene" or not p.genes:
            continue
        names = str(p.genes[0].annotation or p.genes[0].name).split("|")
        hit = next((n for n in landmarks_wanted if n in names), None)
        pl = p.placement(B)
        if hit and pl and not any(l["name"] == hit for l in landmarks):
            near = [i for i, l in enumerate(loci) if l["chrom"] == pl.chrom and abs(l["pos"] - pl.lead_pos) <= 300_000
                    and l["organ"] == landmarks_wanted[hit]]
            landmarks.append({"name": hit, "chrom": pl.chrom, "pos": pl.lead_pos, "organ": landmarks_wanted[hit],
                              "loci_within_300kb": len(near)})

    stats = {"loci": len(loci), "traits": len({(l["panel"], l["trait"]) for l in loci}),
             "eqtls": 44354, "genes": len(symbol), "signals": 480, "records": 49193,
             "known_gene_loci": sum(l["gene"] is not None for l in loci)}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    payload = {"chroms": chroms, "loci": loci, "links": links, "eqtl": eqtl, "bin": BIN, "stats": stats,
               "landmarks": landmarks}
    OUT.write_text("window.ATLAS = " + json.dumps(payload, separators=(",", ":")) + ";\n")
    print(f"{len(loci)} loci, {len(links)} links, {OUT.stat().st_size / 1e3:.0f} kB -> {OUT}")
    print(pd.Series([l["organ"] for l in loci]).value_counts().to_dict())
    print("landmarks:", [(l["name"], l["loci_within_300kb"]) for l in landmarks])


if __name__ == "__main__":
    main()
