# plocust validation (rice, first round: 2026-10-04)

All scripts run from the repository root after `data/download.sh` (data are not committed).
Results are in `results/`.

| | Question | Script | Data |
|---|---|---|---|
| B | Do passports land on the right place in another assembly or build? | `b_anchoring.py` (+ `fetch_orthologs.py`) | funRiceGenes cloned genes; IRGSP-1.0 → MH63, ZS97, N22, Azucena, MSU6 |
| C | Does the same-signal test beat coordinate overlap? | `c_simulation.py` | simulated traits on real 3K genotypes, two panels |
| D | Can loci from two real studies on two builds be matched? | `d_real_gwas.py` | RDP1 (44K, MSU6, published p) vs 3K (core SNPs, IRGSP-1.0, REGENIE) |
| E | Do matched loci recover the known genes? | `d_real_gwas.py` (`genes` step) | rice known-gene DB (4,314 genes) |

## B. Anchoring across assemblies

Gene passports were built on IRGSP-1.0 and placed on each assembly by their three 201 bp anchors.
Truth: Ensembl Compara one-to-one orthologs (300 random cloned genes), or MSU6 gene models via MSU IDs.
"Same coordinates" reuses the IRGSP-1.0 interval unchanged, as a coordinate-only database would.

| Target | Genes | Placed | Correct | Correct when placed | Same coordinates correct |
|---|---|---|---|---|---|
| MH63RS2 (indica) | 256 | 98.8% | 98.4% | 99.6% | 1.2% |
| ZS97RS3 (indica) | 256 | 98.4% | 98.4% | 100% | 2.3% |
| N22RS2 (aus) | 261 | 99.6% | 99.2% | 99.6% | 1.1% |
| AzucenaRS1 (trop. japonica) | 265 | 98.9% | 98.9% | 100% | 1.9% |
| MSU6 (older Nipponbare build) | 272 | 98.9% | 98.9% | 100% | 63.2% |

The 2 wrong placements (of 1,310) were both flagged `split` (anchors disagreed), so neither went undetected.
Across the whole database, 98–99% of the 4,314 genes were placed on each assembly.

**Real old-study anchoring (from D):** 45/45 RDP1 lead SNPs anchored from MSU6 onto IRGSP-1.0
land on exactly the position published in the RDP1 MSU7 SNP table. The MSU6 → IRGSP-1.0 shift of those SNPs is
median 1.1 kb, up to 185 kb.

## C. Simulation: same signal vs nearby distinct signal

3K accessions with phenotypes split by PC1 into panel A (1,409) and panel B (856). Each of 20 replicates has
8 shared causal SNPs, 4 *nearby-distinct* pairs (100–400 kb apart, r² < 0.1, one SNP active in each panel),
and 2 private SNPs per panel; each explains 4% of variance. GWAS: REGENIE + 5 PCs per panel.
All pairs of loci (A × B) with leads within 1 Mb were scored.

| Method | Called same | Precision | Recall | Nearby-distinct pairs wrongly called same |
|---|---|---|---|---|
| **plocust same-signal (r² ≥ 0.5 in panel B)** | 28 | **0.75** | **0.88** | **0 / 12** |
| LD-block interval overlap | 36 | 0.56 | 0.83 | 5 / 12 |
| lead distance ≤ 100 kb | 16 | 0.69 | 0.46 | 0 / 12 |
| lead distance ≤ 250 kb | 26 | 0.58 | 0.63 | 4 / 12 |
| lead distance ≤ 500 kb | 41 | 0.46 | 0.79 | 10 / 12 |

24 truly-same pairs in total. plocust was the only method with both high recall and no false merges of
nearby-distinct signals. Coordinate rules trade one for the other.

Caveats:
- **Small numbers.** Panel B has little power: on average it detects 1.05 of the 8 shared causal SNPs.
  More replicates or a larger second panel are needed before quoting these as final.
- **Precision is conservative.** All 7 "false" same calls have r² 0.53–1.0 between the leads, but one locus
  in each pair was left unassigned by the truth rule (its lead had r² < 0.3 with every causal SNP).
- 157 of 472 detected loci were unassigned; many are extra loci from clump fragmentation in long-range LD.

## D. Real data: RDP1 (MSU6) vs 3K (IRGSP-1.0)

RDP1: published mixed-model p-values, p < 1e-5. 3K: REGENIE (λ 0.91–1.32), Bonferroni.
RDP1 loci were anchored onto IRGSP-1.0, then compared with the same-signal test using 3K LD.
Panel length and flag leaf length were dropped: no RDP1 SNP reaches p < 1e-5.

| Trait (RDP1 / 3K) | RDP1 loci | 3K loci | Best match per RDP1 locus |
|---|---|---|---|
| Plant height / culm length | 2 | 14 | 1 same (identical lead SNP, r² = 1.0), 1 ambiguous |
| Amylose / endosperm type | 6 | 16 | *Wx* locus: same (leads 13 kb apart, r² = 1.0) |
| Pericarp colour / seed coat | 36 | 9 | *Rc* locus: **ambiguous** (leads 0.9 kb apart, r² = 0.13) |
| Awn presence | 1 | 0 | no 3K signal |

**The *Rc* case is the main finding of the first round.** Both studies hit *Rc*, but the two lead SNPs are
nearly uncorrelated in the 3K panel. *Rc* has several functional alleles, and the classic *rc* allele is a
14 bp deletion that is in neither SNP set, so the two panels tag it with different SNPs. Lead-to-lead
r² is therefore not enough at allelic series. plocust says "ambiguous" here, not "distinct", which is the
honest answer, but a haplotype- or credible-set-based test is needed (see open items).

## E. Known genes

| Trait | Expected gene | RDP1 locus | 3K locus |
|---|---|---|---|
| Pericarp colour | *Rc* | in locus (4.9 kb from lead) | in locus (0.2 kb) |
| Amylose / endosperm | *Wx* | in locus (lead inside gene) | in locus (13 kb) |
| Plant height | *sd1* | 270 kb from nearest lead, outside | inside the 0.95 Mb locus, 270 kb from lead |
| Awn | *An-1*, *LABA1*, *GAD1*, ... | not found | not found |

Loci intervals in rice can be wide (3K plant height: 156 cloned genes inside its loci), so the locus card
ranks known genes by **trait-keyword match first**, then distance. *sd1* ranks second for the 3K plant-height
locus although 270 kb from the lead.

Example cards: `results/d_real_gwas/card_3k_*.png`.

## Open items found by this round

1. **Allelic series (*Rc*)**: add a multi-SNP same-signal test (haplotype r², or overlap of anchored credible sets).
2. **Credible sets of anchored passports** are only used on the source build; anchoring the credible-set
   SNPs would let the test use them across builds.
3. **Simulation power**: more replicates and a larger second panel.
4. **Clump fragmentation** in long-range LD inflates locus counts (36 RDP1 pericarp loci); consider
   merging clumps whose leads are in LD or whose LD blocks overlap.
