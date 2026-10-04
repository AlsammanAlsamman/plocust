"""Match passports: across GWAS studies and against the known-loci database.

All comparisons happen on one genome build, so passports from another build
are anchored first (see `anchor`). When a genotype panel on that build is
given, the same-signal test uses LD measured in that panel:

    r2 >= same_r2      -> "same"             (the two leads tag one signal)
    r2 <  distinct_r2  -> "distinct_nearby"  (independent signals in one region)
    otherwise          -> "ambiguous"

Without genotypes the call falls back to distance ("same_by_position").
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .io import Genotypes
from .ld import r2_matrix
from .passport import LocusKind, LocusPassport


@dataclass
class MatchConfig:
    max_distance_kb: float = 1000  # pairs further apart are not compared
    same_r2: float = 0.5
    distinct_r2: float = 0.1
    proxy_kb: float = 10  # a lead missing from the panel is replaced by the nearest SNP within this distance
    position_same_kb: float = 100  # without genotypes: leads this close are called the same
    cs_top: int = 10  # credible-set variants per locus used in the LD test
    gene_pad_kb: float = 50  # a gene within the locus interval +- this is "gene_in_locus"


def _lead(p: LocusPassport, build: str):
    pl = p.placement(build)
    return None if pl is None or pl.lead_pos is None else pl


def _panel_index(geno: Genotypes, chrom: str, positions, proxy_bp: int) -> list[int]:
    out = []
    for pos in positions:
        v = geno.nearest(chrom, int(pos), proxy_bp)
        if v is not None:
            out.append(int(v["idx"]))
    return out


def signal_r2(a: LocusPassport, b: LocusPassport, build: str, geno: Genotypes, cfg: MatchConfig) -> Optional[float]:
    """Highest r2 between the two loci's lead (and, on their own build, credible-set) variants."""
    pa, pb = _lead(a, build), _lead(b, build)

    def positions(p: LocusPassport, pl) -> list[int]:
        pos = [pl.lead_pos]
        if pl.method == "source" and p.credible_set:  # credible-set positions are only valid on the source build
            top = sorted(p.credible_set.variants, key=lambda v: -v.pip)[: cfg.cs_top]
            pos += [v.pos for v in top]
        return list(dict.fromkeys(pos))

    proxy = int(cfg.proxy_kb * 1000)
    ia = _panel_index(geno, pa.chrom, positions(a, pa), proxy)
    ib = _panel_index(geno, pb.chrom, positions(b, pb), proxy)
    if not ia or not ib:
        return None
    g = geno.read(ia + ib)
    return float(np.nanmax(r2_matrix(g[:, : len(ia)], g[:, len(ia) :])))


def match_pair(a: LocusPassport, b: LocusPassport, build: str, geno: Optional[Genotypes] = None,
               cfg: Optional[MatchConfig] = None) -> Optional[dict]:
    """Compare two passports on `build`; None when they are not on it or too far apart."""
    cfg = cfg or MatchConfig()
    pa, pb = _lead(a, build), _lead(b, build)
    if pa is None or pb is None or pa.chrom != pb.chrom:
        return None
    dist = abs(pa.lead_pos - pb.lead_pos)
    if dist > cfg.max_distance_kb * 1000:
        return None
    overlap = pa.start <= pb.end and pb.start <= pa.end
    row = {"query": a.passport_id, "target": b.passport_id, "target_kind": b.kind.value, "chrom": pa.chrom,
           "query_lead": pa.lead_pos, "target_lead": pb.lead_pos, "distance_bp": dist, "overlap": overlap,
           "r2": None, "call": None, "score": None}

    if b.kind is LocusKind.gene or a.kind is LocusKind.gene:
        pad = cfg.gene_pad_kb * 1000
        inside = pa.start - pad <= pb.lead_pos <= pa.end + pad or pb.start - pad <= pa.lead_pos <= pb.end + pad
        row["call"] = "gene_in_locus" if inside else "gene_nearby"
        row["score"] = round(1 - dist / (cfg.max_distance_kb * 1000), 4)
        return row

    r2 = signal_r2(a, b, build, geno, cfg) if geno is not None else None
    if r2 is not None:
        row["r2"] = round(r2, 4)
        row["score"] = round(r2, 4)
        row["call"] = "same" if r2 >= cfg.same_r2 else "distinct_nearby" if r2 < cfg.distinct_r2 else "ambiguous"
    else:
        row["score"] = round(1 - dist / (cfg.max_distance_kb * 1000), 4)
        row["call"] = "same_by_position" if dist <= cfg.position_same_kb * 1000 or overlap else "nearby"
    return row


def compare(query: list[LocusPassport], target: list[LocusPassport], build: str, geno: Optional[Genotypes] = None,
            cfg: Optional[MatchConfig] = None) -> pd.DataFrame:
    """All query x target pairs within max_distance_kb on `build`, with calls."""
    cfg = cfg or MatchConfig()
    by_chrom: dict[str, list] = {}
    for t in target:
        pl = _lead(t, build)
        if pl is not None:
            by_chrom.setdefault(pl.chrom, []).append(t)
    rows = []
    for q in query:
        pl = _lead(q, build)
        if pl is None:
            continue
        for t in by_chrom.get(pl.chrom, []):
            if t.passport_id == q.passport_id and t is q:
                continue
            r = match_pair(q, t, build, geno, cfg)
            if r is not None:
                rows.append(r)
    cols = ["query", "target", "target_kind", "chrom", "query_lead", "target_lead", "distance_bp", "overlap", "r2",
            "call", "score"]
    return pd.DataFrame(rows, columns=cols)


def best_matches(pairs: pd.DataFrame) -> pd.DataFrame:
    """For each query, its best-scoring target among non-gene matches."""
    loci = pairs[pairs["target_kind"] != "gene"]
    if loci.empty:
        return loci
    return loci.sort_values(["query", "score"], ascending=[True, False]).groupby("query", as_index=False).head(1)


def coordinate_overlap(a: LocusPassport, b: LocusPassport, build: str, pad_kb: float = 0) -> bool:
    """Baseline: do the two loci's intervals overlap (after padding) on `build`?"""
    pa, pb = a.placement(build), b.placement(build)
    if pa is None or pb is None or pa.chrom != pb.chrom:
        return False
    pad = pad_kb * 1000
    return pa.start - pad <= pb.end + pad and pb.start - pad <= pa.end + pad
