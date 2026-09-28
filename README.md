# GLP-1 Receptor Agonists and Suicide-Attempt Risk in Youth: A Triangulation Study

Analysis code and derived result tables for the manuscript submitted to the *Journal of the American Academy of Child & Adolescent Psychiatry* (New Research; under review).

## Study design

Three evidence layers with largely non-overlapping bias structures are triangulated:

1. **FAERS pharmacovigilance** — pediatric (<18 y) vs adult disproportionality (ROR, Bayesian IC025) for suicidality and other psychiatric event groups across GLP-1 receptor agonists, 2012Q4–2026Q2 (55 quarters).
2. **YRBS 2025 adolescent epidemiology** — survey-weighted logistic regression of weight-loss intention/behavior on suicidal ideation and attempts in US high-school students.
3. **Developmental drug-target Mendelian randomization** — GLP1R cis-instruments for childhood BMI z-score (MoBa 2022 waves: birth–8 y) and adult BMI (GIANT × Million Veteran Program meta-analysis, N ≈ 1.11 M) against suicide-attempt GWAS (iPSYCH; Models 1–3).

## Repository layout

```
scripts/    analysis and plotting pipeline (Python 3; R/TwoSampleMR replication noted in text)
results/    derived summary tables (CSV) produced by the pipeline
```

Key scripts:

| Script | Purpose |
|---|---|
| `faers_download.py`, `faers_process.py`, `faers_pediatric.py`, `faers_analyze.py` | FAERS quarterly download, deduplication (latest version per CASEID), pediatric subsetting, ROR/IC025 |
| `yrbs_pilot.py`, `yrbs_weightcontrol.py` | YRBS 2025 survey-weighted models |
| `extract_cis.py`, `clump_harmonize.py` | GLP1R cis-window extraction (locus ±100 kb, r² < 0.3) and harmonization |
| `run_child_mr.py`, `run_child_mr_developmental.py` | wave-specific childhood cis-MR and fixed-infant-instrument developmental scan |
| `extract_adult_meta_mvp.py`, `run_adult_meta_mvp_mr.py`, `run_adult_relaxed.py` | adult BMI meta-analysis instruments and MR |
| `run_ld_gls.py` | LD-corrected generalized-least-squares sensitivity |
| `plot_submission_fig1.py`, `plot_developmental_pub.py`, `plot_submission_figS1.py` | manuscript figures (vector PDF) |

## Data sources

- **FAERS**: FDA public quarterly data files (<https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-public-dashboard>)
- **YRBS 2025**: CDC public-use file (<https://www.cdc.gov/yrbs/>)
- **MoBa childhood BMI GWAS summary statistics**: Helgeland et al., *Nature Metabolism* 2022 (Norwegian Institute of Public Health)
- **iPSYCH suicide-attempt GWAS**: available to qualified researchers on application to the iPSYCH consortium (Danish institutional affiliation and Danish Data Protection Agency approval required)
- **GIANT / Million Veteran Program BMI**: public repositories

Raw GWAS summary statistics and FAERS raw quarterly files are not redistributed here; all derived quantities needed to reproduce the reported estimates are in `results/`.

## Status

Manuscript under review. This repository will be made public upon acceptance, with a citable DOI via Zenodo.
