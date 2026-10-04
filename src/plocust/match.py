"""Match passports: across GWAS studies and against the known-loci database.

All comparisons happen on one genome build, so passports from another build
are anchored first (see `anchor`). When a genotype panel on that build is
given, the same-signal test uses LD measured in that panel:

    r2 >= same_r2      -> "same"             (the two leads tag one signal)
    r2 <  distinct_r2  -> "distinct_nearby"  (independent signals in one region)
    otherwise          -> "ambiguous"

When both studies report effects, the direction test asks whether the
trait-increasing alleles sit on the same haplotype in the panel (sign of the
lead-to-lead correlation x the two effect signs). A "same" call with opposite
directions is downgraded to "ambiguous".

Without genotypes the call falls back to distance ("same_by_position").
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .ids import reverse_complement

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
    coloc: bool = True  # also run summary-statistics colocalization when imprints and genotypes allow


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


def _panel_sign(p: LocusPassport, pl, v) -> Optional[int]:
    """Sign of p's effect expressed for the panel's counted allele (A1), or None if unknown."""
    s = p.signal
    if s is None or s.beta is None or s.beta == 0 or not np.isfinite(s.beta):
        return None
    ea = (s.effect_allele or s.lead.alt or "").upper()
    if pl.strand == "-":  # passport alleles are on the other strand of this build
        ea = reverse_complement(ea) if ea and set(ea) <= set("ACGTN") else ea
    sign = 1 if s.beta > 0 else -1
    if ea == str(v["a1"]).upper():
        return sign
    if ea == str(v["a2"]).upper():
        return -sign
    return None


def direction(a: LocusPassport, b: LocusPassport, build: str, geno: Genotypes, min_r2: float = 0.3) -> Optional[str]:
    """'concordant' / 'discordant' effect direction of two leads in LD, or None when it cannot be told."""
    pa, pb = _lead(a, build), _lead(b, build)
    va, vb = geno.nearest(pa.chrom, pa.lead_pos, 0), geno.nearest(pb.chrom, pb.lead_pos, 0)
    if va is None or vb is None:
        return None
    sa, sb = _panel_sign(a, pa, va), _panel_sign(b, pb, vb)
    if sa is None or sb is None:
        return None
    g = geno.read([int(va["idx"]), int(vb["idx"])])
    ok = ~np.isnan(g).any(axis=1)
    if ok.sum() < 10 or g[ok, 0].std() == 0 or g[ok, 1].std() == 0:
        return None
    r = float(np.corrcoef(g[ok, 0], g[ok, 1])[0, 1])
    if r * r < min_r2:
        return None
    return "concordant" if sa * sb * np.sign(r) > 0 else "discordant"


def profile_similarity(a: LocusPassport, b: LocusPassport, build: str) -> Optional[float]:
    """Correlation of the two imprints' significance profiles (oriented to `build`), or None."""
    ia, ib = a.imprint, b.imprint
    if ia is None or ib is None or ia.bins != ib.bins or ia.window_bp != ib.window_bp:
        return None
    pa, pb = np.array(ia.profile), np.array(ib.profile)
    if a.placement(build) and a.placement(build).strand == "-":
        pa = pa[::-1]
    if b.placement(build) and b.placement(build).strand == "-":
        pb = pb[::-1]
    if pa.std() == 0 or pb.std() == 0:
        return None
    return round(float(np.corrcoef(pa, pb)[0, 1]), 4)


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
           "r2": None, "direction": None, "profile_corr": profile_similarity(a, b, build), "call": None,
           "score": None, "PP.H3": None, "PP.H4": None, "coloc_snps": None, "coloc_call": None}

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
        row["direction"] = direction(a, b, build, geno)
        if row["call"] == "same" and row["direction"] == "discordant":
            row["call"] = "ambiguous"
    if geno is not None and cfg.coloc:
        from .coloc import coloc_pair

        res = coloc_pair(a, b, build, geno)
        if res is not None:
            row.update({k: res[k] for k in ("PP.H3", "PP.H4", "coloc_snps", "coloc_call")})
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
            "direction", "profile_corr", "call", "score", "PP.H3", "PP.H4", "coloc_snps", "coloc_call"]
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
