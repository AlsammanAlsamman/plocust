# LOCUST (`plocust`)

**LOCUs Signature Tracker**: build-independent *locus passports* for plant GWAS.

A GWAS locus is usually stored as coordinates on one genome build. Those break
when the build, the marker panel or the population changes. A plocust passport
describes a locus by what it is:

- **sequence anchors** (lead SNP and block edges), which place it on any assembly by alignment
- **the statistical signal** (lead, effect, credible set)
- **LD / haplotypes** in the source panel
- **nearby genes** and **neighbouring loci**

The passport ID comes from the lead-SNP anchor sequence (GA4GH VRS-style
digest, strand-independent), so it stays the same across builds.

**Status:** early development (schema v0.1). First crop: rice. See [PLAN.md](PLAN.md).

## Install (development)

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest
```

## CLI

```bash
plocust schema                 # print the passport JSON Schema
plocust validate p1.json ...   # check passport files
```

Planned: `identify`, `anchor`, `match`, `compare`, `db`.
