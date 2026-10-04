"""Passports from published hit tables: lead SNPs only, no full summary statistics, no genotypes.

Many papers publish only their significant associations (lead SNP, p-value, sometimes effect and
target gene). Each row becomes a passport with a lead anchor from the reference genome, its signal and,
for eQTL tables, the regulated gene. Credible sets, LD blocks and imprints need full summary statistics
and are left empty; matching uses a reference LD panel.
"""

from __future__ import annotations

from typing import Optional

import pandas as pd

from .identify import IdentifyConfig, _check_lead_allele
from .io import Genome, normalize_chrom
from .passport import Anchor, AnchorRole, Gene, LocusKind, LocusPassport, Placement, Signal, Source, Trait, Variant


def passports_from_hits(hits: pd.DataFrame, genome: Genome, species: str, study: str, kind: LocusKind = LocusKind.gwas,
                        genes: Optional[pd.DataFrame] = None, n_samples: Optional[int] = None,
                        cfg: Optional[IdentifyConfig] = None, trait: Optional[str] = None) -> list[LocusPassport]:
    """One passport per row.

    `hits` needs chrom, pos, p; optional: id, ref, alt, beta, trait, gene (target gene, for eQTL).
    `genes` (chrom, start, end, id, name) places the target gene and lists the nearest genes.
    """
    cfg = cfg or IdentifyConfig()
    gene_at = genes.set_index("id") if genes is not None else None
    out = []
    for r in hits.itertuples(index=False):
        chrom, pos = normalize_chrom(r.chrom), int(r.pos)
        seq, off = genome.flank(chrom, pos, cfg.flank)
        if len(seq) < 2 * cfg.flank // 3:
            continue
        ref = getattr(r, "ref", None) if isinstance(getattr(r, "ref", None), str) else seq[off]
        alt = getattr(r, "alt", None) if isinstance(getattr(r, "alt", None), str) else "N"
        anchor = Anchor(role=AnchorRole.lead_snp, sequence=seq, variant_offset=off, label=str(getattr(r, "id", "")) or None)
        notes = ["published hit: no credible set, LD block or imprint"]
        if (msg := _check_lead_allele(anchor, ref, alt)) is not None:
            notes.append(msg)
        gene_list = []
        target = getattr(r, "gene", None)
        if isinstance(target, str) and gene_at is not None and target in gene_at.index:
            g = gene_at.loc[target]
            g = g.iloc[0] if isinstance(g, pd.DataFrame) else g
            dist = 0 if g["start"] <= pos <= g["end"] else int(min(abs(g["start"] - pos), abs(g["end"] - pos)))
            local = g["chrom"] == chrom
            gene_list.append(Gene(id=target, name=g.get("name") if isinstance(g.get("name"), str) else None,
                                  start=int(g["start"]), end=int(g["end"]), distance_bp=dist if local else 10**9,
                                  annotation="target gene" + ("" if local else f" (chr{g['chrom']}, distant)")))
        trait_name = getattr(r, "trait", None) if isinstance(getattr(r, "trait", None), str) else trait
        if kind is LocusKind.eqtl and isinstance(target, str):
            trait_name = f"expression of {target}"
        beta = getattr(r, "beta", None)
        out.append(LocusPassport(
            kind=kind,
            species=species,
            trait=Trait(name=trait_name or "unknown"),
            source=Source(study=study, n_samples=n_samples),
            anchors=[anchor],
            placements=[Placement(build=genome.build, chrom=chrom, start=pos, end=pos, lead_pos=pos)],
            signal=Signal(lead=Variant(id=str(getattr(r, "id", "")) or None, chrom=chrom, pos=pos, ref=ref, alt=alt),
                          beta=float(beta) if beta is not None and pd.notna(beta) else None,
                          pvalue=max(float(r.p), 1e-300)),
            genes=gene_list,
            notes="; ".join(notes),
        ))
    return out
