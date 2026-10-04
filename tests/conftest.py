"""Synthetic genome, second assembly, genotype panel and GWAS for end-to-end tests."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest
from scipy import stats

from plocust.ids import reverse_complement

CHROM_LEN = 300_000
SNP_STEP = 1_000
BLOCK_SNPS = 30
CAUSAL_POS = 100_500
INSERT_AT, INSERT_LEN = 50_000, 5_000


def write_bed(prefix: Path, g: np.ndarray, bim: pd.DataFrame, iids) -> None:
    """Minimal PLINK 1 .bed writer (g: samples x variants, counts of A1, NaN = missing)."""
    n, m = g.shape
    code = np.where(np.isnan(g), 1, np.select([g == 2, g == 1], [0, 2], 3)).astype(np.uint8)
    pad = (-n) % 4
    code = np.vstack([code, np.zeros((pad, m), np.uint8)]).T.reshape(m, -1, 4)
    packed = (code[..., 0] | code[..., 1] << 2 | code[..., 2] << 4 | code[..., 3] << 6).astype(np.uint8)
    with open(f"{prefix}.bed", "wb") as fh:
        fh.write(b"\x6c\x1b\x01" + packed.tobytes())
    bim[["chrom", "id", "cm", "pos", "a1", "a2"]].to_csv(f"{prefix}.bim", sep="\t", header=False, index=False)
    pd.DataFrame({"f": iids, "i": iids, "p": 0, "m": 0, "s": 0, "y": -9}).to_csv(
        f"{prefix}.fam", sep=" ", header=False, index=False)


def write_fasta(path: Path, seqs: dict) -> Path:
    with open(path, "w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return path


@dataclass
class World:
    genome: Path
    target: Path
    target_inv: Path
    gff: Path
    bfile: Path
    sumstats: pd.DataFrame
    sumstats_path: Path
    chrom1: str


@pytest.fixture(scope="session")
def world(tmp_path_factory) -> World:
    d = tmp_path_factory.mktemp("world")
    rng = np.random.default_rng(1)
    bases = np.array(list("ACGT"))
    chr1 = "".join(rng.choice(bases, CHROM_LEN))
    chr2 = "".join(rng.choice(bases, CHROM_LEN))
    genome = write_fasta(d / "src.fa", {"chr01": chr1, "chr02": chr2})
    insert = "".join(rng.choice(bases, INSERT_LEN))
    t1 = chr1[: INSERT_AT - 1] + insert + chr1[INSERT_AT - 1 :]
    target = write_fasta(d / "tgt.fa", {"1": t1, "2": chr2})
    target_inv = write_fasta(d / "tgt_inv.fa", {"1": reverse_complement(chr1), "2": chr2})

    genes = [("chr01", 99_000, 103_000, "+", "G1"), ("chr01", 160_000, 165_000, "-", "G2"),
             ("chr01", 240_000, 242_000, "+", "G3")]
    gff = d / "genes.gff3"
    gff.write_text("##gff-version 3\n" + "".join(
        f"{c}\ttest\tgene\t{s}\t{e}\t.\t{st}\t.\tID=gene:{i};biotype=protein_coding;gene_id={i};Name={i}\n"
        for c, s, e, st, i in genes))

    # Inbred panel: blocks of 30 SNPs; each sample carries one of 4 founder haplotypes per block.
    n = 300
    pos = np.arange(SNP_STEP // 2, CHROM_LEN, SNP_STEP)
    m = len(pos)
    g = np.empty((n, m))
    for b0 in range(0, m, BLOCK_SNPS):
        b1 = min(m, b0 + BLOCK_SNPS)
        founders = rng.integers(0, 2, size=(4, b1 - b0)) * 2
        g[:, b0:b1] = founders[rng.integers(0, 4, n)]
    ref = np.array([chr1[p - 1] for p in pos])
    alt = np.array([{"A": "G", "G": "A", "C": "T", "T": "C"}[b] for b in ref])
    bim = pd.DataFrame({"chrom": "1", "id": [f"s{p}" for p in pos], "cm": 0, "pos": pos, "a1": alt, "a2": ref})
    iids = [f"i{k}" for k in range(n)]
    bfile = d / "panel"
    write_bed(bfile, g, bim, iids)

    causal = int(np.flatnonzero(pos == CAUSAL_POS)[0])
    y = 0.8 * g[:, causal] / 2 + rng.normal(0, 1, n)
    rows = []
    for j in range(m):
        res = stats.linregress(g[:, j], y) if g[:, j].std() > 0 else None
        beta, se, p = (res.slope, res.stderr, res.pvalue) if res else (0.0, np.nan, 1.0)
        rows.append(("chr01", pos[j], f"s{pos[j]}", ref[j], alt[j], alt[j], beta, se, max(p, 1e-300), n))
    ss = pd.DataFrame(rows, columns=["CHR", "BP", "SNP", "REF", "ALT", "A1", "BETA", "SE", "P", "N"])
    ss_path = d / "gwas.tsv"
    ss.to_csv(ss_path, sep="\t", index=False)
    return World(genome, target, target_inv, gff, bfile, ss, ss_path, "1")
