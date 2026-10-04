"""Readers for GWAS summary statistics, genomes, gene annotation and genotypes."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .ids import normalize_sequence

# ---------------------------------------------------------------- chromosomes

_CHROM_PREFIX = re.compile(r"^(chr|chromosome|chrom)_?", re.IGNORECASE)


def normalize_chrom(chrom) -> str:
    """'chr01', 'Chr1', 'chromosome1', 'chr01|13101', 1 -> '1'. Other names are kept as they are."""
    s = _CHROM_PREFIX.sub("", re.split(r"[|\s]", str(chrom).strip())[0])
    return str(int(s)) if s.isdigit() else s


# ---------------------------------------------------------- summary statistics

SUMSTATS_COLUMNS = ["chrom", "pos", "id", "ref", "alt", "effect_allele", "beta", "se", "p", "af", "n"]

# Lower-case column name -> standard name, for the common GWAS tools.
_ALIASES = {
    "chrom": ["chrom", "#chrom", "chr", "chromosome", "chr_id"],
    "pos": ["pos", "bp", "position", "ps", "base_pair_location", "genpos"],
    "id": ["id", "snp", "rs", "marker", "snp_id", "variant_id", "rsid"],
    "ref": ["ref", "allele0", "a2", "other_allele"],
    "alt": ["alt", "allele1", "a1", "effect_allele"],
    "effect_allele": ["a1", "allele1", "effect_allele"],
    "beta": ["beta", "effect", "b"],
    "se": ["se", "stderr", "standard_error"],
    "p": ["p", "pvalue", "p_value", "p_wald", "p_lrt", "p.value", "pval", "p_score"],
    "af": ["af", "maf", "a1_freq", "a1freq", "freq", "effect_allele_frequency"],
    "n": ["n", "obs_ct", "n_miss_complement", "nobs"],
}


def read_sumstats(path: str | Path, columns: Optional[dict[str, str]] = None, sep: Optional[str] = None) -> pd.DataFrame:
    """Read GWAS results into the standard columns (see SUMSTATS_COLUMNS).

    Column names of PLINK2 --glm, GEMMA, GAPIT and most CSV/TSV outputs are
    recognised. `columns` maps standard names to the file's names when they
    are not (e.g. {"p": "MLM_pvalue"}).
    """
    raw = pd.read_csv(path, sep=sep, engine="python" if sep is None else "c", comment=None)
    return standardize_sumstats(raw, columns)


def standardize_sumstats(raw: pd.DataFrame, columns: Optional[dict[str, str]] = None) -> pd.DataFrame:
    lower = {c.lower().strip(): c for c in raw.columns}
    found: dict[str, str] = {}
    for std, aliases in _ALIASES.items():
        if columns and std in columns:
            found[std] = columns[std]
            continue
        hit = next((lower[a] for a in aliases if a in lower), None)
        if hit is not None:
            found[std] = hit
    if "p" not in found and "log10p" in lower:  # REGENIE and others report -log10(p)
        raw = raw.assign(_p=np.maximum(10.0 ** -pd.to_numeric(raw[lower["log10p"]], errors="coerce"), 1e-300))
        found["p"] = "_p"
    missing = [c for c in ("chrom", "pos", "p") if c not in found]
    if missing:
        raise ValueError(f"summary statistics need columns {missing}; pass `columns` to map them")

    out = pd.DataFrame({std: raw[col].values for std, col in found.items()})
    for c in SUMSTATS_COLUMNS:
        if c not in out:
            out[c] = np.nan
    out["chrom"] = out["chrom"].map(normalize_chrom)
    out["pos"] = out["pos"].astype(int)
    out["p"] = pd.to_numeric(out["p"], errors="coerce")
    out = out[out["p"].notna() & (out["p"] > 0) & (out["p"] <= 1)]
    if out["id"].isna().all():
        out["id"] = out["chrom"] + ":" + out["pos"].astype(str)
    out["id"] = out["id"].astype(str)
    return out[SUMSTATS_COLUMNS].sort_values(["chrom", "pos"]).reset_index(drop=True)


# ------------------------------------------------------------------- genomes


class Genome:
    """Indexed FASTA (plain or bgzip) with tolerant chromosome names."""

    def __init__(self, path: str | Path, build: Optional[str] = None):
        import pysam

        self.path = Path(path)
        self.build = build or self.path.name.split(".dna")[0]
        self._fa = pysam.FastaFile(str(self.path))
        self._names = {normalize_chrom(n): n for n in self._fa.references}

    def name(self, chrom) -> str:
        key = normalize_chrom(chrom)
        if key not in self._names:
            raise KeyError(f"chromosome {chrom!r} not in {self.path.name}")
        return self._names[key]

    def length(self, chrom) -> int:
        return self._fa.get_reference_length(self.name(chrom))

    def fetch(self, chrom, start: int, end: int) -> str:
        """1-based inclusive interval, clipped to the chromosome."""
        name = self.name(chrom)
        start = max(1, start)
        end = min(self._fa.get_reference_length(name), end)
        if end < start:  # interval entirely past the chromosome end
            return ""
        return self._fa.fetch(name, start - 1, end).upper()

    def flank(self, chrom, pos: int, flank: int) -> tuple[str, int]:
        """Sequence of pos +- flank, and the 0-based offset of pos inside it."""
        start = max(1, pos - flank)
        return self.fetch(chrom, start, pos + flank), pos - start


class FlankTable:
    """Per-SNP flanking sequences, for studies whose genome build is not at hand.

    Expects columns for the SNP ID and the sequence; the variant is written
    as [A/G] (Illumina style) or as a single IUPAC code in the middle.
    """

    _BRACKET = re.compile(r"\[([ACGTN-]+)/([ACGTN-]+)\]", re.IGNORECASE)

    def __init__(self, table: pd.DataFrame, id_col: str, seq_col: str):
        self._seqs = dict(zip(table[id_col].astype(str), table[seq_col].astype(str)))

    @classmethod
    def read(cls, path: str | Path, id_col: Optional[str] = None, seq_col: Optional[str] = None, sep: str = "\t"):
        t = pd.read_csv(path, sep=sep)
        id_col = id_col or t.columns[0]
        seq_col = seq_col or next(c for c in t.columns if "seq" in c.lower() or "flank" in c.lower())
        return cls(t, id_col, seq_col)

    def __contains__(self, snp_id: str) -> bool:
        return str(snp_id) in self._seqs

    def flank(self, snp_id: str, ref: Optional[str] = None) -> tuple[str, int]:
        """Flank with the variant replaced by `ref` (or its first allele); returns (sequence, offset)."""
        raw = self._seqs[str(snp_id)]
        m = self._BRACKET.search(raw)
        if m:
            allele = (ref or m.group(1)).upper()
            left, right = raw[: m.start()], raw[m.end() :]
            seq = left + allele + right
            return normalize_sequence(seq.replace("-", "")), len(left)
        # IUPAC code at the variant: take the first ambiguous base
        for i, b in enumerate(raw.upper()):
            if b not in "ACGTN":
                allele = (ref or "N").upper()
                return normalize_sequence(raw[:i] + allele + raw[i + 1 :]), i
        raise ValueError(f"no variant marked in flank of {snp_id}")


# -------------------------------------------------------------------- genes

_GFF_ATTR = re.compile(r"([^=;]+)=([^;]*)")


def read_genes(path: str | Path, feature: str = "gene") -> pd.DataFrame:
    """Genes from a GFF3 (plain or gzip): chrom, start, end, strand, id, name, biotype."""
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            attrs = dict(_GFF_ATTR.findall(f[8]))
            gid = attrs.get("gene_id") or attrs.get("ID", "").split(":")[-1]
            rows.append((normalize_chrom(f[0]), int(f[3]), int(f[4]), f[6], gid, attrs.get("Name"), attrs.get("biotype")))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "id", "name", "biotype"])


# ---------------------------------------------------------------- genotypes

_BED_LOOKUP = np.array(
    # PLINK .bed 2-bit codes: 00 hom A1, 01 missing, 10 het, 11 hom A2 -> count of A1
    [[2, np.nan, 1, 0][(byte >> (2 * k)) & 3] for byte in range(256) for k in range(4)],
    dtype=np.float32,
).reshape(256, 4)


@dataclass
class Genotypes:
    """PLINK 1 binary fileset (.bed/.bim/.fam), read lazily by region."""

    prefix: Path
    variants: pd.DataFrame  # chrom, id, pos, a1, a2, idx
    samples: pd.DataFrame  # fid, iid
    _keep: Optional[np.ndarray] = None

    @classmethod
    def open(cls, prefix: str | Path) -> Genotypes:
        prefix = Path(str(prefix).removesuffix(".bed"))
        bim = pd.read_csv(f"{prefix}.bim", sep=r"\s+", header=None, names=["chrom", "id", "cm", "pos", "a1", "a2"],
                          dtype={"chrom": str, "id": str, "a1": str, "a2": str})
        bim["chrom"] = bim["chrom"].map(normalize_chrom)
        bim["idx"] = np.arange(len(bim))
        fam = pd.read_csv(f"{prefix}.fam", sep=r"\s+", header=None, usecols=[0, 1], names=["fid", "iid"], dtype=str)
        with open(f"{prefix}.bed", "rb") as fh:
            if fh.read(3) != b"\x6c\x1b\x01":
                raise ValueError(f"{prefix}.bed is not a SNP-major PLINK bed file")
        return cls(prefix, bim.drop(columns="cm"), fam)

    @property
    def n_samples(self) -> int:
        return len(self.samples) if self._keep is None else int(self._keep.sum())

    def subset(self, iids) -> Genotypes:
        keep = self.samples["iid"].isin(set(map(str, iids))).to_numpy()
        return Genotypes(self.prefix, self.variants, self.samples, keep)

    def region(self, chrom, start: int, end: int) -> pd.DataFrame:
        chrom = normalize_chrom(chrom)
        v = self.variants
        return v[(v["chrom"] == chrom) & (v["pos"] >= start) & (v["pos"] <= end)]

    def read(self, idx) -> np.ndarray:
        """Matrix samples x variants of A1 counts (NaN = missing)."""
        idx = np.asarray(idx, dtype=np.int64)
        n = len(self.samples)
        nbytes = (n + 3) // 4
        mm = np.memmap(f"{self.prefix}.bed", dtype=np.uint8, mode="r", offset=3)
        packed = mm.reshape(-1, nbytes)[idx]  # variants x bytes
        g = _BED_LOOKUP[packed].reshape(len(idx), nbytes * 4)[:, :n].T
        return g[self._keep] if self._keep is not None else g

    def nearest(self, chrom, pos: int, max_dist: int) -> Optional[pd.Series]:
        r = self.region(chrom, pos - max_dist, pos + max_dist)
        if r.empty:
            return None
        return r.iloc[int(np.argmin(np.abs(r["pos"].to_numpy() - pos)))]
