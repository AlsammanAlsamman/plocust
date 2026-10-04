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
| 7 | `07_finemap_loci.py` | fine-mapping of every GWAS locus (SuSiE, in-sample LD adjusted for the GWAS PCs) |
| 8 | `08_cross_subpop.py` | cross-subpopulation fine-mapping on 12 real 3K traits, dense SNPs |
| 9 | `09_known_genes_4k.py` | known cloned genes vs fine-mapped credible sets, RiceVarMap 4K genotypes incl. deletions |

## The database (v0.2.0-atlas, 112 MB, local; every GWAS record carries its fine-mapping)

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
| 7 | Fine-mapping with raw LD on mixed-model z-scores invented signals (1,146 "signals"; 88% of regions with LD outliers) | LD adjusted for the GWAS PCs + signals must reach p < 1e-5: 480 signals, 63% regions with outliers |
| 8 | ...but in the 529-line panel the PC adjustment hid *GS3* (PIP 0.21 vs 0.92 with raw LD) | both LD versions are reported; choosing LD per panel is an open problem |
| 9 | Variants sharing a position (SNP + indel) broke lookups | lookups by position no longer assume uniqueness |

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

## Fine-mapping every GWAS locus (step 7)

358 regions, 525 clumped loci → **480 fine-mapped signals** (LD adjusted for the GWAS PCs; a signal needs p < 1e-5).
Most traits now have *fewer* signals than clumps (glutinous endosperm 16 → 6, seed coat 9 → 7), as expected when
one signal had been split. Blade pubescence (263 → 257) stays implausibly high and is flagged.
Caveat: the 3K **core** set is LD-pruned, so its credible sets look artificially small (median 1 SNP); variant-level
resolution needs dense SNPs (steps 8–9).

## Cross-subpopulation fine-mapping on real traits (step 8): dense 3K SNPs, 87 regions, 12 traits

| Method | Signals | Median top credible set | Top set hits the trait's cloned gene (25 regions with one) |
|---|---|---|---|
| full panel (n≈2,100) | 62 | 2 | 7 |
| A alone | 34 | 2 | 8 |
| B alone | 10 | 3 | 4 |
| A + B joint | 40 | 2 | 8 |

On real traits the joint model gains little over the full panel (8 vs 7 genes): the large gain seen in simulation
(step 6, credible sets 90 → 47) does not carry over here, partly because the full panel already contains both
subpopulations. Not yet a reportable improvement.

## Known cloned genes (step 9): 529 accessions, 4K genotypes with deletions, 16 gene × trait tests with a signal

| Method | Credible set hits the gene | Gene variant with PIP ≥ 0.5 |
|---|---|---|
| GWAS lead inside the gene (baseline) | 2/16 | n/a |
| fine-mapping, raw in-sample LD | **3/16** | **2/16** |
| fine-mapping, PC-adjusted LD | 1/16 | 0/16 |

- ***sd1***: found by every method.
- ***GS3***: the gene carries a signal as strong as the window lead (4e-12 vs 3e-12). With raw LD the model puts
  **PIP 0.92 on chr3:16,733,441 (G/T) inside *GS3***, consistent with the known nonsense mutation (to be confirmed
  against the published position); PC-adjusted LD spreads it over 2,209 variants.
- *GW5*, *GW2*, *GL7*, *GW8*: strong signals nearby, but weak inside the genes in this panel: other loci, or causal
  variants (large deletions) that the imputed panel does not represent well.

**Honest summary:** fine-mapping is useful for counting signals, but its variant-level resolution on real rice data is
still limited and depends on how LD is specified. That, not the matching, is now the main scientific open problem.

## Next

1. LD specification per panel: kinship-aware LD (match the mixed model) instead of raw vs PC-adjusted.
2. Joint fine-mapping across *independent* studies (e.g. 3K + 529 panel), where the cross-subpopulation gain should appear.
3. Structural variants (large deletions) in the reference panel, for *GW5*-type causal variants.
4. Tissue-matched expression data for a fair eQTL benchmark.
