# Rice locus atlas: running plocust on public rice data

Scripts run in order from the repository root (data are fetched into `data/`, which is not committed).
Small result tables are in [`results/`](results/).

| Step | Script | What it does |
|---|---|---|
| 1 | `01_gwas.py` | mixed-model GWAS (REGENIE) for 64 traits: 29 from 3K, 35 from RDP1 |
| 2 | `02_eqtl_passports.py` | 44,354 published leaf eQTL hits (Liu et al. 2022) → passports (kind `eqtl`) |
| 3 | `03_gwas_passports.py` | passports for every GWAS trait; RDP1 anchored from MSU6 to IRGSP-1.0 |
| 4 | `04_database.py` | rice passport database: genes + GWAS + eQTL records, chip markers, LD panel |
| 5 | `05_learn.py` | learning analyses: fragmentation, known genes, cross-study matching, GWAS × eQTL |
| 6 | `06_finemap_sim.py` | fine-mapping model on real genotypes: single vs cross-subpopulation, borrowed LD |

## The database (v0.1.0-atlas, 112 MB, local)

| Records | Count |
|---|---|
| cloned genes (funRiceGenes) | 4,314 |
| GWAS loci (3K: 415, RDP1: 110) | 525 |
| leaf eQTLs (Liu et al. 2022) | 44,354 |
| **total records** | **49,193** |
| SNP-chip markers (RDP1 44K) | 36,775 |

A *record* is one locus in one study for one trait. Records of the same DNA share a `passport_id`, so
pleiotropy is a lookup: 2,885 SNPs regulate more than one gene; 25 GWAS lead SNPs act on more than one trait
(e.g. apiculus and leaf-sheath colour, both anthocyanin).

## What the runs taught us, and what changed in the tool

| # | Finding | Change |
|---|---|---|
| 1 | REGENIE over-corrects in the small RDP1 panel (λ 0.37–0.72) | RDP1 uses the published EMMA p-values; REGENIE only for effect signs |
| 2 | One signal split into many loci by long-range LD (blade pubescence 323, culm length 17) | merge clumps whose leads have r² ≥ 0.5 (3K loci 501 → 415); the rest needs fine-mapping |
| 3 | Comparisons too slow on a 4.8M-SNP panel (did not finish in 40 min) | per-chromosome binary-search index, vectorised SNP matching, proxy cache: whole analysis 91 s |
| 4 | One SNP with several traits/genes collided in the database (3,814 records lost) | records keyed by locus + study + trait; `passport_id` stays the locus identity |
| 5 | The LD call was overwritten when colocalization was switched off | bug fixed, regression test added |
| 6 | Published hit tables (no summary statistics) were not importable | `plocust.hits`: passports from lead SNP + p (+ target gene) |

## Fine-mapping model (step 6): real 3K genotypes, 40 regions, 68 causal SNPs

| Method | Power | Correct credible sets | Median credible-set size |
|---|---|---|---|
| subpopulation A alone (n≈1,400) | 0.63 | 0.89 | 90 SNPs |
| subpopulation B alone (n≈850) | 0.35 | 0.85 | 64.5 SNPs |
| **A + B jointly (cross-subpopulation)** | **0.69** | **0.87** | **47 SNPs** |
| A's z-scores with B's LD (borrowed) | 0.06 | 0.02 | 171 false signals |

- **Cross-subpopulation fine-mapping halves credible sets** while keeping them correct: the main methodological result so far.
- **Borrowed LD from another subpopulation is disastrous.** For studies without genotypes, the LD reference must
  match the study population, and the consistency check (`finemap.ld_consistency`) must pass first.

## Cross-study matching (3K vs RDP1)

| Pair | Result |
|---|---|
| plant height ~ culm length (*sd1* region) | same signal |
| amylose ~ glutinous endosperm (*Wx*) | same signal, **opposite effect** (opposite trait coding), colocalized |
| pericarp ~ seed coat colour (*Rc*) | same signal (2 loci), colocalized |
| leaf ~ blade pubescence | 1 same, 3 ambiguous: needs fine-mapping |

## GWAS × leaf eQTL: inconclusive with this data

163 of 525 GWAS loci share a signal with at least one local leaf eQTL. At the 38 loci with a trait-matched cloned gene,
the eQTL gene was never the cloned gene (nearest gene: 5/38). But **only 3 of those 18 cloned genes have any leaf
eQTL**: pigment, grain and pubescence genes act in other tissues. The benchmark cannot be scored with leaf data;
it needs tissue-matched eQTL (panicle, seed) or leaf-acting traits. The 125 "new candidate genes" are not yet trustworthy.

## Next

1. Apply the fine-mapping model to the GWAS loci (replaces pairwise clump merging; resolves pubescence).
2. RiceVarMap 4K genotypes (with deletions) as a second LD panel and for the 533-panel traits (heading, yield, grain size).
3. Known causal variants (*Wx*, *sd1*, *GS3*, *Rc*, *GW5*, *Hd1*, *Ghd7*): do credible sets contain them?
4. Tissue-matched expression data for a fair eQTL benchmark.
