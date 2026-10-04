"""Which genome build do these coordinates belong to?

A breeder often has positions and alleles but not the genome file, or not the
build name. On the right build, the reference base at almost every SNP is one of
its two alleles. On a wrong build the base at that position is unrelated to the
SNP, so the match falls to chance (about 2 bases out of 4, i.e. ~50%).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd

from .io import Genome

_COMP = str.maketrans("ACGT", "TGCA")


@dataclass
class BuildMatch:
    build: str
    tested: int
    match: float  # share of SNPs whose reference base is one of the alleles
    match_flipped: float  # same, allowing the alleles to be on the other strand
    missing_chrom: int

    @property
    def best(self) -> float:
        return max(self.match, self.match_flipped)


def allele_match(ss: pd.DataFrame, genome: Genome, n: int = 2000, seed: int = 0) -> BuildMatch:
    """Share of (biallelic, single-base) SNPs whose reference base is one of their alleles."""
    snv = ss[ss["ref"].astype(str).str.len().eq(1) & ss["alt"].astype(str).str.len().eq(1)
             & ss["ref"].astype(str).str.upper().isin(list("ACGT")) & ss["alt"].astype(str).str.upper().isin(list("ACGT"))]
    # A/T and C/G SNPs match either strand, so they say nothing about the strand
    pair = (snv["ref"].str.upper() + snv["alt"].str.upper())
    snv = snv[~pair.isin(["AT", "TA", "CG", "GC"])]
    if len(snv) > n:
        snv = snv.sample(n, random_state=seed)
    hit = flip = missing = tested = 0
    for chrom, pos, a, b in zip(snv["chrom"], snv["pos"], snv["ref"].str.upper(), snv["alt"].str.upper()):
        try:
            base = genome.fetch(chrom, int(pos), int(pos))
        except KeyError:
            missing += 1
            continue
        if not base:  # past the chromosome end: counts as a mismatch
            tested += 1
            continue
        tested += 1
        hit += base in (a, b)
        flip += base in (a.translate(_COMP), b.translate(_COMP))
    return BuildMatch(genome.build, tested, round(hit / max(tested, 1), 4), round(flip / max(tested, 1), 4), missing)


def detect_build(ss: pd.DataFrame, genomes: Mapping[str, Genome] | list[Genome], n: int = 2000) -> pd.DataFrame:
    """Rank candidate builds by allele match; the top row is the likely build."""
    gs = list(genomes.values()) if isinstance(genomes, Mapping) else list(genomes)
    rows = []
    for g in gs:
        m = allele_match(ss, g, n=n)
        rows.append({"build": m.build, "snps_tested": m.tested, "allele_match": m.match,
                     "allele_match_other_strand": m.match_flipped, "chrom_missing": m.missing_chrom})
    df = pd.DataFrame(rows)
    df["score"] = np.maximum(df["allele_match"], df["allele_match_other_strand"])
    df = df.sort_values("score", ascending=False).reset_index(drop=True)
    df["verdict"] = ""
    if len(df):
        top = df.loc[0, "score"]
        second = df.loc[1, "score"] if len(df) > 1 else 0.0
        if top >= 0.97 and top - second >= 0.05:
            df.loc[0, "verdict"] = "likely build"
        elif top >= 0.97:
            df.loc[0, "verdict"] = "fits, but so do others"
        else:
            df.loc[0, "verdict"] = "no candidate fits well"
    return df
