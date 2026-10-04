"""Place passports on any genome assembly by aligning their sequence anchors."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .io import normalize_chrom
from .passport import AnchorRole, LocusPassport, Placement, PlacementFlag


class Aligner:
    """minimap2 (via mappy) index of an assembly. The index is cached as <fasta>.sr.mmi."""

    def __init__(self, fasta: str | Path, build: Optional[str] = None, preset: str = "sr", threads: int = 4):
        import mappy

        fasta = Path(fasta)
        mmi = fasta.with_name(fasta.name + f".{preset}.mmi")
        if mmi.exists():
            self._a = mappy.Aligner(str(mmi), preset=preset, best_n=5)
        else:
            self._a = mappy.Aligner(str(fasta), preset=preset, best_n=5, n_threads=threads, fn_idx_out=str(mmi))
        if not self._a:
            raise RuntimeError(f"could not build a minimap2 index for {fasta}")
        self.build = build or fasta.name.split(".dna")[0]

    def map(self, seq: str):
        return list(self._a.map(seq))


@dataclass
class AnchorHit:
    chrom: str
    pos: int  # 1-based position of the anchor's variant_offset on the target
    strand: int  # +1 / -1 relative to the anchor as stored
    mapq: int
    identity: float


def _project(hit, offset: int, qlen: int) -> int:
    """0-based target position aligned to query position `offset`, walking the CIGAR."""
    if hit.strand == 1:
        q, r, target = hit.q_st, hit.r_st, offset
    else:  # alignment is of the reverse complement
        q, r, target = qlen - hit.q_en, hit.r_st, qlen - 1 - offset
    if target < q:  # inside a clipped start: extrapolate
        return r - (q - target)
    for length, op in hit.cigar:
        if op in (0, 7, 8):  # M, =, X
            if q + length > target:
                return r + (target - q)
            q, r = q + length, r + length
        elif op == 1:  # insertion: query only
            if q + length > target:
                return r
            q += length
        elif op in (2, 3):  # deletion / skip: target only
            r += length
    return r + (target - q)


def map_anchor(aligner: Aligner, seq: str, offset: int, min_mapq: int = 20) -> tuple[Optional[AnchorHit], int]:
    """Best unique hit of one anchor, and the number of hits found."""
    hits = aligner.map(seq)
    if not hits:
        return None, 0
    best = next((h for h in hits if h.is_primary), hits[0])
    if best.mapq < min_mapq:
        return None, len(hits)
    hit = AnchorHit(
        chrom=normalize_chrom(best.ctg),
        pos=_project(best, offset, len(seq)) + 1,
        strand=best.strand,
        mapq=best.mapq,
        identity=round(best.mlen / max(best.blen, 1), 4),
    )
    return hit, len(hits)


def place(passport: LocusPassport, aligner: Aligner, min_mapq: int = 20, max_span_factor: float = 3.0,
          min_span_slack: int = 200_000) -> Optional[Placement]:
    """Placement of a passport on the aligner's assembly, or None when the primary anchor cannot be placed."""
    hits: dict[AnchorRole, Optional[AnchorHit]] = {}
    multi = False
    for a in passport.anchors:
        hit, n = map_anchor(aligner, a.sequence, a.variant_offset, min_mapq)
        hits[a.role] = hit
        multi |= hit is None and n > 0

    primary = hits.get(passport.primary_role)
    if primary is None:
        return None

    flags: set[PlacementFlag] = set()
    if any(h is None for h in hits.values()):
        flags.add(PlacementFlag.missing_anchor)
    if multi:
        flags.add(PlacementFlag.multi_mapping)
    if primary.strand == -1:
        flags.add(PlacementFlag.inverted)

    src = next((p for p in passport.placements if p.method == "source"), None)
    src_span = (src.end - src.start) if src else 0
    max_span = max(min_span_slack, max_span_factor * src_span)
    agree = {}
    for role, h in hits.items():
        if h is None:
            continue
        if h.chrom != primary.chrom or abs(h.pos - primary.pos) > max_span:
            flags.add(PlacementFlag.split)
        else:
            agree[role] = h

    left, right = agree.get(AnchorRole.block_left), agree.get(AnchorRole.block_right)
    if left and right:
        # On the + strand the left edge should come first; on the - strand the order flips.
        if (left.pos > right.pos) != (primary.strand == -1):
            flags.add(PlacementFlag.out_of_order)
    positions = [h.pos for h in agree.values()]
    return Placement(
        build=aligner.build,
        chrom=primary.chrom,
        start=max(1, min(positions)),
        end=max(positions),
        strand="+" if primary.strand == 1 else "-",
        method="anchor",
        lead_pos=primary.pos,
        anchors_placed=len(agree),
        anchors_total=len(passport.anchors),
        flags=sorted(flags, key=lambda f: f.value),
    )


def anchor_passports(passports: list[LocusPassport], aligner: Aligner, **kw) -> list[Optional[Placement]]:
    """Place each passport and add the placement to it (replacing any earlier one for that build)."""
    out = []
    for p in passports:
        pl = place(p, aligner, **kw)
        p.placements = [x for x in p.placements if x.build != aligner.build or x.method == "source"]
        if pl is not None and p.placement(aligner.build) is None:
            p.placements.append(pl)
        out.append(pl)
    return out


def check_uniqueness(passports: list[LocusPassport], aligner: Aligner, min_mapq: int = 20) -> None:
    """Set Anchor.unique by mapping each anchor back to its own genome."""
    for p in passports:
        for a in p.anchors:
            hit, _ = map_anchor(aligner, a.sequence, a.variant_offset, min_mapq)
            a.unique = hit is not None
