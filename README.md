# LOCUST (`plocust`)

**LOCUs Signature Tracker**: build-independent *locus passports* for plant GWAS.

A GWAS locus is usually stored as coordinates on one genome build. Those break
when the build, the marker panel or the population changes. A plocust passport
describes a locus by what it is:

- **sequence anchors** (lead SNP and LD-block edges, 201 bp each) that place it on any assembly by alignment
- **the statistical signal** (lead, effect, p, credible set)
- **LD and haplotypes** in the source panel
- **nearby genes** and **neighbouring loci**

The passport ID (`OsLP.<digest>`) is computed from the lead-SNP anchor sequence
(GA4GH VRS-style digest, strand-independent), so it never changes with the build.

**Status:** early development (schema 0.1.0). First crop: rice. See [PLAN.md](PLAN.md)
and the [validation results](validation/README.md).

## Install

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest
```

Requires Python ≥ 3.10. `mappy` (minimap2) and `pysam` are installed from PyPI.

## Workflow

```bash
# 1. GWAS results -> passports (PLINK2, GEMMA, GAPIT, REGENIE or any CSV/TSV)
plocust identify --sumstats gwas.tsv --genome IRGSP-1.0.fa --build IRGSP-1.0 \
    --gff IRGSP-1.0.gff3.gz --bfile panel --species "Oryza sativa" \
    --trait "plant height" --ontology TO:0000207 --study my_gwas --check-unique --out loci.jsonl

#    old study without its genome at hand? give per-SNP flanking sequences instead:
plocust identify --sumstats old.tsv --flanks snp_flanks.tsv --build MSU6 ... --out old.jsonl

# 2. place passports on another assembly or build, by sequence
plocust anchor old.jsonl --genome IRGSP-1.0.fa --build IRGSP-1.0 --out old.irgsp.jsonl

# 3. compare two studies (same-signal test uses LD in the panel you give)
plocust compare old.irgsp.jsonl loci.jsonl --build IRGSP-1.0 --bfile panel --out pairs.tsv

# 4. match to the known-loci database
plocust match loci.jsonl --db plocust-db-rice.sqlite --build IRGSP-1.0 --out known.tsv

# 5. one figure per locus for breeders
plocust card loci.jsonl --index 0 --sumstats gwas.tsv --gff IRGSP-1.0.gff3.gz \
    --bfile panel --db plocust-db-rice.sqlite --out locus.png
```

`identify` writes a `.jsonl` (one passport per line) and a `.tsv` summary next to it.

### Calls made by `compare` / `match`

| Call | Meaning |
|---|---|
| `same` | leads (or credible-set SNPs) in LD, r² ≥ 0.5 in the panel given |
| `distinct_nearby` | r² < 0.1: independent signals in the same region |
| `ambiguous` | 0.1 ≤ r² < 0.5 |
| `same_by_position` / `nearby` | no genotypes given: decided by distance only |
| `gene_in_locus` / `gene_nearby` | a known gene inside the locus interval (±50 kb), or within range |

Placement flags from `anchor`: `missing_anchor`, `multi_mapping`, `split`
(anchors disagree), `inverted`, `out_of_order`.

## The known-loci database

The database is a separate, versioned SQLite file. Releases will be fetched with
`plocust db download rice`; none is published yet. Build it locally from
funRiceGenes:

```bash
plocust db build-genes --genes geneInfo.table.txt --keywords geneKeyword.table.txt \
    --gff IRGSP-1.0.gff3.gz --genome IRGSP-1.0.fa --build IRGSP-1.0 \
    --assembly MH63RS2.fa:MH63RS2 --assembly ZS97RS3.fa:ZS97RS3 \
    --species "Oryza sativa" --version 0.1.0 --out plocust-db-rice.sqlite
plocust db info plocust-db-rice.sqlite
plocust db query plocust-db-rice.sqlite --build IRGSP-1.0 --trait "plant height"
```

## Python API

```python
from plocust.io import read_sumstats, Genome, Genotypes, read_genes
from plocust.identify import identify_loci, StudyInfo
from plocust.anchor import Aligner, anchor_passports
from plocust.match import compare

loci = identify_loci(read_sumstats("gwas.tsv"), StudyInfo("Oryza sativa", "plant height", "my_gwas", "IRGSP-1.0"),
                     genome=Genome("IRGSP-1.0.fa"), geno=Genotypes.open("panel"))
anchor_passports(loci, Aligner("MH63RS2.fa", build="MH63RS2"))
```
