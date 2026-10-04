"""F2 at scale: 300 random RDP1 chip markers, anchored onto IRGSP-1.0 from manifest flanks only."""

from pathlib import Path

import pandas as pd

from plocust.anchor import Aligner, map_anchor
from plocust.io import FlankTable

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "validation" / "results" / "f_new_features"

info = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU6.SNP_Information.MSU7.txt", sep="\t", index_col=False)
truth = dict(zip(info["SNPID"], zip(info["CHR"].astype(str), pd.to_numeric(info["position.IRGSP1-MSU7"], errors="coerce"))))
f7 = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU7.SNP_flanking_seq.txt", sep="\t", dtype=str)
f6 = pd.read_csv(DATA / "rdp1/RiceDiversity.44K.MSU6.SNP_flanking_seq.txt", sep="\t", dtype=str)
ids = f7["snp_id"][f7["snp_id"].isin(f6["snp_id"]) & f7["snp_id"].isin(truth)].sample(300, random_state=1)
aligner = Aligner(DATA / "ref/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa", build="IRGSP-1.0")
rows = []
for label, t, a, b in (("41 bp", f7, "msu7_seq1", "msu7_seq2"), ("33 bp", f6, "X5p_MSU6", "X3p_MSU6")):
    t = t.assign(sequence=t[a] + "[" + t["alleles"] + "]" + t[b])
    flanks = FlankTable(t, "snp_id", "sequence")
    for sid in ids:
        seq, off = flanks.flank(sid)
        hit, n = map_anchor(aligner, seq, off)
        chrom, pos = truth[sid]
        rows.append({"flank": label, "snp": sid, "placed": hit is not None,
                     "exact": hit is not None and hit.chrom == chrom and hit.pos == pos, "hits": n})
df = pd.DataFrame(rows)
df.to_csv(OUT / "f2_markers_300.tsv", sep="\t", index=False)
s = df.groupby("flank").agg(markers=("snp", "size"), placed=("placed", "mean"), exact=("exact", "mean"))
s["wrong_when_placed"] = df[df["placed"]].groupby("flank")["exact"].apply(lambda x: int((~x).sum()))
print(s.round(3).to_string())
