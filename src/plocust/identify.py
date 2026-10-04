"""Turn GWAS summary statistics into locus passports."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import norm

from .io import FlankTable, Genome, Genotypes
from .ld import r2_matrix, r2_with
from .passport import (
    Anchor,
    AnchorRole,
    CredibleSet,
    CredibleVariant,
    Gene,
    Imprint,
    LDBlock,
    LocusPassport,
    Neighbour,
    Placement,
    Signal,
    Source,
    Trait,
    Variant,
)


@dataclass
class IdentifyConfig:
    p_threshold: Optional[float] = None  # None = Bonferroni 0.05 / number of tests
    window_kb: float = 500  # clumping radius around a lead
    clump_r2: float = 0.1  # with genotypes: SNPs with r2 >= this join the lead's clump
    merge_r2: Optional[float] = 0.5  # merge clumps whose leads are in LD (one signal split by long-range LD); None = off
    merge_kb: float = 2000  # only clumps this close are considered for merging
    min_significant: int = 1  # drop clumps with fewer significant SNPs (singletons are often artefacts)
    flank: int = 100  # anchor = variant +- flank bp
    block_r2: float = 0.5  # LD block = region around the lead where r2 with it stays >= this
    block_max_kb: float = 1500
    block_max_gap_kb: float = 50  # the block ends after this distance without a SNP at r2 >= block_r2
    haplotype_snps: int = 30
    haplotype_min_freq: float = 0.02
    inversion_min_kb: float = 300
    inversion_min_r2: float = 0.9
    credible_coverage: float = 0.95
    abf_prior_sd: float = 0.2  # prior effect sd, in phenotype sd units (Wakefield)
    abf_default_ratio: float = 0.95  # W/(V+W) when n and af are not available
    gene_window_kb: float = 50
    max_genes: int = 15
    neighbour_kb: float = 2000
    imprint_kb: float = 250  # imprint window: lead +- this
    imprint_bins: int = 25
    imprint_nodes: int = 40  # top significant SNPs in the LD network
    imprint_edge_r2: float = 0.5


@dataclass
class StudyInfo:
    species: str
    trait: str
    study: str
    build: str
    trait_ontology: Optional[str] = None
    panel: Optional[str] = None
    n_samples: Optional[int] = None
    gwas_method: Optional[str] = None
    doi: Optional[str] = None
    extra: dict = field(default_factory=dict)


# ------------------------------------------------------------------ statistics


def z_scores(ss: pd.DataFrame) -> np.ndarray:
    """|z| from the two-sided p-value, the same quantity that ranks the leads.

    beta/se is used only where p underflowed (<= 1e-300); tests such as Firth or LRT give p-values
    that differ from the Wald ratio, and mixing the two would let the credible set disagree with the lead.
    """
    p = ss["p"].to_numpy(dtype=float)
    z = np.abs(norm.isf(p / 2))
    if "beta" in ss and "se" in ss:
        bz = np.abs(ss["beta"].to_numpy(dtype=float) / ss["se"].to_numpy(dtype=float))
        under = (p <= 1e-300) & np.isfinite(bz)
        z[under] = np.maximum(z[under], bz[under])
    return z


def abf_pips(ss: pd.DataFrame, cfg: IdentifyConfig) -> np.ndarray:
    """Posterior inclusion probabilities under one causal variant (Wakefield approximate Bayes factors)."""
    z = z_scores(ss)
    n = ss["n"].to_numpy(dtype=float) if "n" in ss else np.full(len(ss), np.nan)
    af = ss["af"].to_numpy(dtype=float) if "af" in ss else np.full(len(ss), np.nan)
    maf = np.minimum(af, 1 - af)
    v = 1.0 / (2 * n * maf * (1 - maf))  # sampling variance of beta in phenotype-sd units
    w = cfg.abf_prior_sd**2
    r = np.where(np.isfinite(v) & (v > 0), w / (v + w), cfg.abf_default_ratio)
    log_abf = 0.5 * np.log1p(-r) + 0.5 * r * z**2
    log_abf -= log_abf.max()
    pip = np.exp(log_abf)
    return pip / pip.sum()


def credible_set(ss: pd.DataFrame, cfg: IdentifyConfig) -> CredibleSet:
    pip = abf_pips(ss, cfg)
    order = np.argsort(-pip)
    k = int(np.searchsorted(np.cumsum(pip[order]), cfg.credible_coverage) + 1)
    keep = order[: min(k, len(order))]
    return CredibleSet(
        method="wakefield_abf",
        coverage=cfg.credible_coverage,
        variants=[
            CredibleVariant(id=ss["id"].iat[i], pos=int(ss["pos"].iat[i]), pip=round(float(pip[i]), 6)) for i in keep
        ],
    )


# --------------------------------------------------------------------- clumping


@dataclass
class _Clump:
    lead: int  # row index into the summary statistics
    members: list[int]


def clump(ss: pd.DataFrame, threshold: float, cfg: IdentifyConfig, geno: Optional[Genotypes] = None) -> list[_Clump]:
    """Greedy clumping by p-value: distance only, or distance + LD when genotypes are given."""
    sig = ss.index[ss["p"] <= threshold].to_numpy()
    sig = sig[np.argsort(ss.loc[sig, "p"].to_numpy(), kind="stable")]
    chroms = ss.loc[sig, "chrom"].to_numpy()
    pos = ss.loc[sig, "pos"].to_numpy()
    remaining = np.ones(len(sig), dtype=bool)
    window = cfg.window_kb * 1000
    clumps = []
    for k, lead in enumerate(sig):
        if not remaining[k]:
            continue
        near_k = np.flatnonzero(remaining & (chroms == chroms[k]) & (np.abs(pos - pos[k]) <= window))
        if geno is not None and len(near_k) > 1:
            r2 = _r2_to_lead(ss, lead, sig[near_k], geno)
            near_k = near_k[(sig[near_k] == lead) | ~np.isfinite(r2) | (r2 >= cfg.clump_r2)]
        remaining[near_k] = False
        members = sig[near_k][np.argsort(pos[near_k], kind="stable")]
        clumps.append(_Clump(int(lead), [int(i) for i in members]))
    if geno is not None and cfg.merge_r2 is not None:
        clumps = merge_clumps(ss, clumps, geno, cfg)
    return [c for c in clumps if len(c.members) >= cfg.min_significant]


def merge_clumps(ss: pd.DataFrame, clumps: list[_Clump], geno: Genotypes, cfg: IdentifyConfig) -> list[_Clump]:
    """Merge clumps whose leads are in LD (r2 >= merge_r2): one signal split by long-range LD.

    Clumps are visited strongest first; a clump joins the first stronger clump it is in LD with.
    """
    if len(clumps) < 2:
        return clumps
    leads = [c.lead for c in clumps]
    idx = _variant_index(ss, leads, geno)
    chroms = ss.loc[leads, "chrom"].to_numpy()
    pos = ss.loc[leads, "pos"].to_numpy()
    owner = list(range(len(clumps)))
    for j in range(1, len(clumps)):
        if idx[j] < 0:
            continue
        cand = [i for i in range(j) if owner[i] == i and chroms[i] == chroms[j] and idx[i] >= 0
                and abs(pos[i] - pos[j]) <= cfg.merge_kb * 1000]
        if not cand:
            continue
        g = geno.read([idx[j]] + [idx[i] for i in cand])
        r2 = r2_with(g[:, 1:], g[:, 0])
        best = int(np.argmax(r2))
        if r2[best] >= cfg.merge_r2:
            owner[j] = cand[best]
    merged: dict[int, list[int]] = {}
    for j, o in enumerate(owner):
        merged.setdefault(o, []).extend(clumps[j].members)
    return [_Clump(clumps[o].lead, sorted(set(m), key=lambda i: ss.at[i, "pos"])) for o, m in merged.items()]


def _variant_index(ss: pd.DataFrame, rows, geno: Genotypes) -> np.ndarray:
    """Genotype column index for each summary-statistics row (-1 if not genotyped), matched by chrom+pos."""
    key = geno.variants.set_index(["chrom", "pos"])["idx"]
    key = key[~key.index.duplicated()]
    idx = pd.MultiIndex.from_arrays([ss.loc[rows, "chrom"], ss.loc[rows, "pos"]])
    return key.reindex(idx).fillna(-1).astype(int).to_numpy()


def _r2_to_lead(ss, lead, rows, geno) -> np.ndarray:
    idx = _variant_index(ss, [lead] + list(rows), geno)
    out = np.full(len(rows), np.nan)
    if idx[0] < 0:
        return out
    ok = idx[1:] >= 0
    if ok.any():
        g = geno.read(np.concatenate([[idx[0]], idx[1:][ok]]))
        out[ok] = r2_with(g[:, 1:], g[:, 0])
    return out


# ------------------------------------------------------------------ imprint


def signed_z(ss: pd.DataFrame) -> tuple[np.ndarray, bool]:
    """z for the alt allele (sign from beta) when effects are reported, else |z|."""
    z = np.abs(norm.isf(ss["p"].to_numpy(dtype=float) / 2))
    beta = pd.to_numeric(ss["beta"], errors="coerce").to_numpy() if "beta" in ss else np.full(len(ss), np.nan)
    if not np.isfinite(beta).any():
        return z, False
    sign = np.sign(np.nan_to_num(beta))
    ea = ss["effect_allele"].astype(str).str.upper().to_numpy()
    ref = ss["ref"].astype(str).str.upper().to_numpy()
    sign = np.where(ea == ref, -sign, sign)  # express every effect for the alt allele
    return z * sign, True


def imprint(ss: pd.DataFrame, lead_row, threshold: float, other_leads: list[int], cfg: IdentifyConfig,
            geno: Optional[Genotypes] = None) -> Imprint:
    """`other_leads`: positions of the study's other leads on the same chromosome."""
    chrom, lead_pos = lead_row["chrom"], int(lead_row["pos"])
    w = int(cfg.imprint_kb * 1000)
    r = ss[(ss["chrom"] == chrom) & (ss["pos"] >= lead_pos - w) & (ss["pos"] <= lead_pos + w)]
    z, signed = signed_z(r)
    az = np.abs(z)
    edges = np.linspace(lead_pos - w, lead_pos + w, cfg.imprint_bins + 1)
    which = np.clip(np.searchsorted(edges, r["pos"].to_numpy(), side="right") - 1, 0, cfg.imprint_bins - 1)
    profile = np.zeros(cfg.imprint_bins)
    np.maximum.at(profile, which, az)
    lp = -np.log10(r["p"].to_numpy())
    half = r["pos"].to_numpy()[lp >= lp.max() / 2]
    sig = np.flatnonzero(r["p"].to_numpy() <= threshold)

    ld_edges = []
    if geno is not None and len(sig) > 1:
        top = sig[np.argsort(-az[sig])][: cfg.imprint_nodes]
        idx = _variant_index(r, r.index[top], geno)
        ok = idx >= 0
        if ok.sum() > 1:
            g = geno.read(idx[ok])
            r2 = r2_matrix(g, g)
            nodes = top[ok]
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    if r2[i, j] >= cfg.imprint_edge_r2:
                        ld_edges.append((int(nodes[i]), int(nodes[j]), round(float(r2[i, j]), 3)))
    left = [x for x in other_leads if x < lead_pos]
    right = [x for x in other_leads if x > lead_pos]
    territory = (max(lead_pos - w, (max(left) + lead_pos) // 2 + 1) if left else lead_pos - w,
                 min(lead_pos + w, (min(right) + lead_pos) // 2) if right else lead_pos + w)
    n_loci_in_window = 1 + sum(abs(x - lead_pos) <= w for x in other_leads)
    return Imprint(
        window_bp=w, threshold=threshold, z_signed=signed,
        territory_start=max(1, territory[0]), territory_end=territory[1],
        ids=r["id"].astype(str).tolist(), pos=r["pos"].astype(int).tolist(),
        ref=r["ref"].fillna("N").astype(str).tolist(), alt=r["alt"].fillna("N").astype(str).tolist(),
        z=[round(float(x), 3) for x in z], bins=cfg.imprint_bins, profile=[round(float(x), 3) for x in profile],
        half_max_width_bp=int(half.max() - half.min()) if len(half) else 0, n_significant=len(sig),
        independent_signals=n_loci_in_window, ld_edge_r2=cfg.imprint_edge_r2 if geno is not None else None,
        ld_edges=ld_edges,
    )


# ----------------------------------------------------------- LD block / haplotypes


def ld_block(chrom: str, lead_pos: int, geno: Genotypes, cfg: IdentifyConfig) -> Optional[tuple[LDBlock, pd.DataFrame]]:
    """Extend from the lead while r2 with it stays >= block_r2 (tolerating short dips)."""
    half = int(cfg.block_max_kb * 1000 / 2)
    region = geno.region(chrom, lead_pos - half, lead_pos + half)
    lead_rows = region[region["pos"] == lead_pos]
    if lead_rows.empty:
        return None
    g = geno.read(region["idx"].to_numpy())
    li = region.index.get_loc(lead_rows.index[0])
    r2 = r2_with(g, g[:, li])
    pos = region["pos"].to_numpy()
    gap = cfg.block_max_gap_kb * 1000
    left = right = li
    for i in range(li - 1, -1, -1):
        if pos[left] - pos[i] > gap:
            break
        if r2[i] >= cfg.block_r2:
            left = i
    for i in range(li + 1, len(r2)):
        if pos[i] - pos[right] > gap:
            break
        if r2[i] >= cfg.block_r2:
            right = i

    inside = np.arange(left, right + 1)
    strong = inside[r2[inside] >= cfg.block_r2]
    start, end = int(region["pos"].iat[left]), int(region["pos"].iat[right])
    n_hap, freqs = _haplotypes(g[:, strong], cfg)
    length_kb = (end - start + 1) / 1000
    block = LDBlock(
        start=start,
        end=end,
        method=f"r2>={cfg.block_r2} with lead",
        r2_threshold=cfg.block_r2,
        n_haplotypes=n_hap,
        haplotype_freqs=freqs,
        inversion_like=bool(length_kb >= cfg.inversion_min_kb and np.median(r2[strong]) >= cfg.inversion_min_r2),
    )
    edges = region.iloc[[left, right]]
    return block, edges


def _haplotypes(g: np.ndarray, cfg: IdentifyConfig) -> tuple[Optional[int], list[float]]:
    """Haplotypes over block SNPs, treating homozygous calls as haplotypes (inbred panels)."""
    if g.shape[1] < 2:  # one SNP is an allele, not a haplotype
        return None, []
    cols = np.linspace(0, g.shape[1] - 1, min(cfg.haplotype_snps, g.shape[1])).round().astype(int)
    h = g[:, np.unique(cols)]
    h = h[~np.isnan(h).any(axis=1) & (h != 1).all(axis=1)]
    if len(h) == 0:
        return None, []
    _, counts = np.unique(h, axis=0, return_counts=True)
    f = np.sort(counts / counts.sum())[::-1]
    f = f[f >= cfg.haplotype_min_freq]
    return int(len(f)), [round(float(x), 4) for x in f[:10]]


# ------------------------------------------------------------------ anchors


def _anchor(role: AnchorRole, chrom, pos, snp_id, ref, cfg, genome, flanks) -> Optional[Anchor]:
    if flanks is not None and snp_id is not None and str(snp_id) in flanks:
        seq, off = flanks.flank(snp_id, ref if isinstance(ref, str) else None)
    elif genome is not None:
        seq, off = genome.flank(chrom, int(pos), cfg.flank)
    else:
        return None
    return Anchor(role=role, sequence=seq, variant_offset=off, label=str(snp_id) if snp_id is not None else None)


def _check_lead_allele(anchor: Anchor, ref, alt) -> Optional[str]:
    base = anchor.sequence[anchor.variant_offset]
    alleles = {a.upper() for a in (ref, alt) if isinstance(a, str) and len(a) == 1}
    if alleles and base not in alleles:
        return f"lead anchor base {base} matches neither allele {sorted(alleles)}: wrong genome build?"
    return None


# -------------------------------------------------------------------- genes


def nearby_genes(genes: pd.DataFrame, chrom, lead_pos, start, end, cs: Optional[CredibleSet], cfg) -> list[Gene]:
    pad = int(cfg.gene_window_kb * 1000)
    g = genes[(genes["chrom"] == chrom) & (genes["end"] >= start - pad) & (genes["start"] <= end + pad)].copy()
    if g.empty:
        return []
    g["dist"] = np.maximum(0, np.maximum(g["start"] - lead_pos, lead_pos - g["end"]))
    cs_pos = np.array([v.pos for v in cs.variants]) if cs else np.array([])
    out = []
    for _, r in g.sort_values("dist").head(cfg.max_genes).iterrows():
        in_cs = bool(((cs_pos >= r["start"]) & (cs_pos <= r["end"])).any())
        out.append(Gene(id=r["id"], name=r["name"] if isinstance(r["name"], str) else None, start=int(r["start"]),
                        end=int(r["end"]), distance_bp=int(r["dist"]), in_credible_set=in_cs,
                        annotation=r["biotype"] if isinstance(r["biotype"], str) else None))
    return out


# ------------------------------------------------------------------- driver


def identify_loci(
    ss: pd.DataFrame,
    study: StudyInfo,
    genome: Optional[Genome] = None,
    flanks: Optional[FlankTable] = None,
    genes: Optional[pd.DataFrame] = None,
    geno: Optional[Genotypes] = None,
    cfg: Optional[IdentifyConfig] = None,
) -> list[LocusPassport]:
    """Detect loci in summary statistics and describe each one as a passport.

    Anchors come from `flanks` (per-SNP flanking sequences) when the SNP is
    listed there, otherwise from `genome`. One of the two is required.
    """
    if genome is None and flanks is None:
        raise ValueError("need a genome FASTA or a flanking-sequence table to build anchors")
    cfg = cfg or IdentifyConfig()
    threshold = cfg.p_threshold if cfg.p_threshold is not None else 0.05 / len(ss)
    clumps = clump(ss, threshold, cfg, geno)

    claimed = {i for c in clumps for i in c.members}
    passports = []
    for c in clumps:
        lead = ss.loc[c.lead]
        chrom, lead_pos = lead["chrom"], int(lead["pos"])
        notes = []
        members = ss.loc[c.members]

        block, edges = None, None
        if geno is not None:
            res = ld_block(chrom, lead_pos, geno, cfg)
            if res is None:
                notes.append("lead not in genotype panel; LD block not computed")
            else:
                block, edges = res
        start = min(int(members["pos"].min()), block.start if block else lead_pos)
        end = max(int(members["pos"].max()), block.end if block else lead_pos)

        # Credible set over the tested SNPs in the clump window, minus the members of other clumps
        # and SNPs closer to another lead (one causal variant per clump).
        w = cfg.window_kb * 1000
        win = ss[(ss["chrom"] == chrom) & (ss["pos"] >= lead_pos - w) & (ss["pos"] <= lead_pos + w)]
        win = win[~win.index.isin(claimed.difference(c.members))]
        others = [int(ss.at[o.lead, "pos"]) for o in clumps if o is not c and ss.at[o.lead, "chrom"] == chrom]
        if others:
            d_self = np.abs(win["pos"].to_numpy() - lead_pos)
            d_other = np.min(np.abs(win["pos"].to_numpy()[:, None] - np.array(others)[None, :]), axis=1)
            win = win[d_self <= d_other]
        cs = credible_set(win, cfg)

        anchors = []
        lead_anchor = _anchor(AnchorRole.lead_snp, chrom, lead_pos, lead["id"], lead["ref"], cfg, genome, flanks)
        if lead_anchor is None:
            notes.append(f"no flanking sequence for lead {lead['id']}; locus skipped")
            continue
        if (msg := _check_lead_allele(lead_anchor, lead["ref"], lead["alt"])) is not None:
            notes.append(msg)
        anchors.append(lead_anchor)
        if edges is not None:
            left, right = edges.iloc[0], edges.iloc[1]
        else:
            left, right = members.iloc[0], members.iloc[-1]
        for role, row in ((AnchorRole.block_left, left), (AnchorRole.block_right, right)):
            if int(row["pos"]) == lead_pos:
                continue
            a = _anchor(role, chrom, int(row["pos"]), row["id"], row.get("ref", row.get("a2")), cfg, genome, flanks)
            if a is not None:
                anchors.append(a)

        signal = Signal(
            lead=Variant(id=lead["id"], chrom=chrom, pos=lead_pos,
                         ref=lead["ref"] if isinstance(lead["ref"], str) else "N",
                         alt=lead["alt"] if isinstance(lead["alt"], str) else "N"),
            effect_allele=lead["effect_allele"] if isinstance(lead["effect_allele"], str) else None,
            effect_allele_freq=float(lead["af"]) if pd.notna(lead["af"]) else None,
            beta=float(lead["beta"]) if pd.notna(lead["beta"]) else None,
            se=float(lead["se"]) if pd.notna(lead["se"]) else None,
            pvalue=float(lead["p"]),
            peak_start=int(members["pos"].min()),
            peak_end=int(members["pos"].max()),
            n_significant=len(members),
        )
        passports.append(
            LocusPassport(
                species=study.species,
                trait=Trait(name=study.trait, ontology_id=study.trait_ontology),
                source=Source(study=study.study, panel=study.panel, n_samples=study.n_samples,
                              gwas_method=study.gwas_method, doi=study.doi),
                anchors=anchors,
                placements=[Placement(build=study.build, chrom=chrom, start=start, end=end, lead_pos=lead_pos)],
                signal=signal,
                credible_set=cs,
                ld_block=block,
                imprint=imprint(ss, lead, threshold, others, cfg, geno),
                genes=nearby_genes(genes, chrom, lead_pos, start, end, cs, cfg) if genes is not None else [],
                notes="; ".join(notes) or None,
            )
        )
    add_neighbours(passports, study.build, geno, cfg)
    return passports


def add_neighbours(passports: list[LocusPassport], build: str, geno: Optional[Genotypes], cfg: IdentifyConfig) -> None:
    """Link each locus to other loci of the same study within neighbour_kb."""
    placed = [(p, p.placement(build)) for p in passports]
    placed = [(p, pl) for p, pl in placed if pl is not None and pl.lead_pos is not None]
    for p, pl in placed:
        near = [(q, ql) for q, ql in placed if q is not p and ql.chrom == pl.chrom
                and abs(ql.lead_pos - pl.lead_pos) <= cfg.neighbour_kb * 1000]
        out = []
        for q, ql in near:
            r2 = lead_r2(geno, pl.chrom, pl.lead_pos, ql.lead_pos) if geno is not None else None
            out.append(Neighbour(passport_id=q.passport_id, distance_bp=abs(ql.lead_pos - pl.lead_pos),
                                 r2_with_lead=None if r2 is None else round(r2, 4),
                                 independent=None if r2 is None else r2 < cfg.clump_r2))
        p.neighbours = sorted(out, key=lambda n: n.distance_bp)


def lead_r2(geno: Genotypes, chrom, pos_a: int, pos_b: int, proxy_bp: int = 0) -> Optional[float]:
    """r2 in `geno` between the variants at two positions (or their nearest proxies within proxy_bp)."""
    a = geno.nearest(chrom, pos_a, proxy_bp)
    b = geno.nearest(chrom, pos_b, proxy_bp)
    if a is None or b is None:
        return None
    g = geno.read([int(a["idx"]), int(b["idx"])])
    return float(r2_matrix(g[:, :1], g[:, 1:])[0, 0])
