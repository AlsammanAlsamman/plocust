"""The locus card: one figure that summarises a passport for breeders."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .io import Genotypes
from .ld import r2_with
from .passport import LocusPassport

# Reference palette (light): sequential blue for r2, orange accent for the lead / known genes.
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#a3a29c"
GRID = "#e6e5e0"
SURFACE = "#fcfcfb"
ACCENT = "#eb6834"
R2_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def locus_card(
    passport: LocusPassport,
    out: str | Path,
    sumstats: Optional[pd.DataFrame] = None,
    genes: Optional[pd.DataFrame] = None,
    geno: Optional[Genotypes] = None,
    known: Optional[list[LocusPassport]] = None,
    build: Optional[str] = None,
    p_threshold: Optional[float] = None,
    pad_kb: float = 150,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    src = next((p for p in passport.placements if p.method == "source"), passport.placements[0])
    build = build or src.build
    pl = passport.placement(build) or src
    chrom, lead = pl.chrom, pl.lead_pos
    lo, hi = max(1, pl.start - int(pad_kb * 1000)), pl.end + int(pad_kb * 1000)
    mb = 1e6

    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2, "xtick.color": INK_2,
                         "ytick.color": INK_2, "font.family": "DejaVu Sans"})
    fig = plt.figure(figsize=(13, 6.8), facecolor=SURFACE)
    gs = fig.add_gridspec(3, 2, width_ratios=[2.1, 1], height_ratios=[3, 1.2, 1.1], hspace=0.12, wspace=0.16,
                          left=0.06, right=0.98, top=0.88, bottom=0.08)
    ax = fig.add_subplot(gs[0, 0])
    axg = fig.add_subplot(gs[1, 0], sharex=ax)
    axh = fig.add_subplot(gs[2, 0])
    axi = fig.add_subplot(gs[:, 1])
    for a in (ax, axg, axh, axi):
        a.set_facecolor(SURFACE)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)

    title = f"{passport.trait.name}  ·  {passport.source.study}"
    fig.text(0.06, 0.95, title, fontsize=14, color=INK, weight="bold")
    fig.text(0.06, 0.915, f"{passport.passport_id}   |   chr{chrom}:{pl.start:,}-{pl.end:,} ({build})",
             fontsize=9, color=INK_2, family="DejaVu Sans Mono")

    # --- association panel
    if passport.ld_block and pl.method == "source":
        ax.axvspan(passport.ld_block.start / mb, passport.ld_block.end / mb, color=GRID, zorder=0, lw=0)
    if sumstats is not None and pl.method == "source":
        r = sumstats[(sumstats["chrom"] == chrom) & (sumstats["pos"] >= lo) & (sumstats["pos"] <= hi)]
        y = -np.log10(r["p"].to_numpy())
        x = r["pos"].to_numpy() / mb
        cmap = LinearSegmentedColormap.from_list("r2", R2_RAMP)
        if geno is not None:
            r2 = _r2_to_lead(r, lead, chrom, geno)
            known_r2 = np.isfinite(r2)
            ax.scatter(x[~known_r2], y[~known_r2], s=12, color=MUTED, lw=0, zorder=2)
            sc = ax.scatter(x[known_r2], y[known_r2], c=r2[known_r2], cmap=cmap, vmin=0, vmax=1, s=16, lw=0, zorder=3)
            cax = ax.inset_axes([1.01, 0.05, 0.012, 0.9])  # inset, so the panel keeps the gene track's width
            cb = fig.colorbar(sc, cax=cax)
            cb.set_label("r² with lead", color=INK_2)
            cb.outline.set_visible(False)
        else:
            ax.scatter(x, y, s=14, color=R2_RAMP[4], lw=0, zorder=3)
        if passport.credible_set:
            cs = {v.pos for v in passport.credible_set.variants}
            m = r["pos"].isin(cs).to_numpy()
            ax.scatter(x[m], y[m], s=46, facecolor="none", edgecolor=INK, lw=0.9, zorder=4, label="credible set")
        if p_threshold:
            ax.axhline(-np.log10(p_threshold), color=MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
        if passport.signal:
            ax.scatter([lead / mb], [-np.log10(passport.signal.pvalue)], marker="D", s=70, color=ACCENT,
                       edgecolor=SURFACE, lw=1.5, zorder=5, label="lead SNP")
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, frameon=False, fontsize=8)
    else:
        ax.text(0.5, 0.5, f"placed on {build} by sequence anchors\n(association shown on source build {src.build})",
                ha="center", va="center", color=INK_2, transform=ax.transAxes)
    ax.set_ylabel("−log10 p")
    ax.grid(axis="y", color=GRID, lw=0.6)
    plt.setp(ax.get_xticklabels(), visible=False)

    # --- gene track
    known_ids = {k.genes[0].id: k for k in (known or []) if k.genes}
    ranked = rank_known(passport, list(known_ids.values()), build, lead)
    labelled = {k.genes[0].id for k in ranked[:6]}  # dense regions are unreadable when every gene is labelled
    if genes is not None:
        g = genes[(genes["chrom"] == chrom) & (genes["end"] >= lo) & (genes["start"] <= hi)].sort_values("start")
        rows_end = []
        for _, row in g.iterrows():
            lane = next((i for i, e in enumerate(rows_end) if e < row["start"] - 15_000), len(rows_end))
            if lane == len(rows_end):
                rows_end.append(0)
            rows_end[lane] = row["end"]
            hit = row["id"] in known_ids
            axg.plot([row["start"] / mb, row["end"] / mb], [-lane, -lane], lw=6 if hit else 4,
                     color=ACCENT if hit else MUTED, solid_capstyle="round")
            if row["id"] in labelled:
                sym = known_ids[row["id"]].genes[0].name or row["id"]
                axg.text((row["start"] + row["end"]) / 2 / mb, -lane + 0.45, sym, ha="center", fontsize=8,
                         color=INK, weight="bold")
        axg.set_ylim(-max(len(rows_end), 1) - 0.2, 1.2)
    axg.set_yticks([])
    axg.spines["left"].set_visible(False)
    axg.axvline(lead / mb, color=ACCENT, lw=1, alpha=0.6)
    axg.set_xlabel(f"chr{chrom} position (Mb, {build})")
    ax.set_xlim(lo / mb, hi / mb)

    # --- haplotypes
    blk = passport.ld_block
    if blk and blk.haplotype_freqs:
        f = blk.haplotype_freqs
        axh.barh(range(len(f)), f, color=R2_RAMP[4], height=0.7)
        axh.set_yticks(range(len(f)), [f"H{i + 1}" for i in range(len(f))])
        axh.invert_yaxis()
        for i, v in enumerate(f):
            axh.text(v + 0.01, i, f"{v:.0%}", va="center", fontsize=8, color=INK_2)
        axh.set_xlim(0, 1.1)
        axh.set_xlabel("haplotype frequency in source panel")
        axh.grid(axis="x", color=GRID, lw=0.6)
    else:
        axh.axis("off")

    # --- info panel
    axi.axis("off")
    lines = []
    s = passport.signal
    if s:
        lines += ["SIGNAL", f"lead {s.lead.id}  {s.lead.ref}>{s.lead.alt}", f"p = {s.pvalue:.2e}"]
        if s.beta is not None:
            lines.append(f"effect ({s.effect_allele or s.lead.alt}) = {s.beta:+.3g}")
        if passport.credible_set:
            lines.append(f"credible set: {len(passport.credible_set.variants)} SNPs "
                         f"({passport.credible_set.coverage:.0%})")
    if blk:
        if blk.end > blk.start:
            desc = f"{blk.length_kb:,.1f} kb, {blk.n_haplotypes or '-'} haplotypes"
        else:
            desc = f"lead only (no SNP at r²≥{blk.r2_threshold})"
        lines += ["", "LD BLOCK", desc + ("  [inversion-like]" if blk.inversion_like else "")]
    lines += ["", "ANCHORS"]
    for a in passport.anchors:
        u = {True: "unique", False: "not unique", None: "unchecked"}[a.unique]
        lines.append(f"{a.role.value:<12} {len(a.sequence)} bp  {u}")
    lines += ["", "PLACEMENTS"]
    for p in passport.placements:
        fl = f" [{', '.join(f.value for f in p.flags)}]" if p.flags else ""
        lines.append(f"{p.build:<12} chr{p.chrom}:{(p.lead_pos or p.start):,} {p.strand}{fl}")
    if known:
        lines += ["", "KNOWN GENES (trait match first)"]
        for k in ranked[:6]:
            kp = k.placement(build)
            d = abs(kp.lead_pos - lead) / 1000 if kp else float("nan")
            star = "*" if trait_match(passport, k) else " "
            traits = ", ".join([k.trait.name] + [t.name for t in k.traits_other][:2])
            lines.append(f"{star}{(k.genes[0].name if k.genes else k.passport_id)[:11]:<11} {d:4.0f} kb  {traits[:22]}")
        lines.append("* trait keyword matches this locus")
    axi.text(0.02, 1.0, "\n".join(lines), va="top", ha="left", fontsize=8.2, color=INK, family="DejaVu Sans Mono",
             transform=axi.transAxes, linespacing=1.45)

    out = Path(out)
    fig.savefig(out, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    return out


def trait_match(passport: LocusPassport, known: LocusPassport) -> bool:
    """Does any word of the locus trait (e.g. 'plant height') appear in the known locus's traits?"""
    words = {w for w in passport.trait.name.lower().replace("/", " ").split() if len(w) > 3}
    names = " ".join([known.trait.name] + [t.name for t in known.traits_other]).lower()
    return bool(words) and any(w in names for w in words)


def rank_known(passport: LocusPassport, known: list[LocusPassport], build: str, lead: int) -> list[LocusPassport]:
    """Known loci ordered by trait match, then distance to the lead."""
    def key(k):
        kp = k.placement(build)
        return (not trait_match(passport, k), abs(kp.lead_pos - lead) if kp else float("inf"))
    return sorted(known, key=key)


def _r2_to_lead(r: pd.DataFrame, lead: int, chrom: str, geno: Genotypes) -> np.ndarray:
    v = geno.region(chrom, int(r["pos"].min()), int(r["pos"].max()))
    v = v[~v["pos"].duplicated()]
    idx = v.set_index("pos")["idx"].reindex(r["pos"].to_numpy())
    out = np.full(len(r), np.nan)
    if lead not in v["pos"].to_numpy():
        return out
    lead_idx = int(v.loc[v["pos"] == lead, "idx"].iat[0])
    ok = idx.notna().to_numpy()
    if ok.any():
        g = geno.read(np.concatenate([[lead_idx], idx[ok].astype(int).to_numpy()]))
        out[ok] = r2_with(g[:, 1:], g[:, 0])
    return out
