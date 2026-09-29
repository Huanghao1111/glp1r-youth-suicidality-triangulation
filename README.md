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
```

Derived result tables and manuscript figures are archived on Zenodo: <https://doi.org/10.5281/zenodo.23008872>

Key scripts:

| Script | Purpose |
|---|---|
| `faers_download.py`, `faers_process.py`, `faers_pediatric.py`, `faers_analyze.py` | FAERS quarterly download, deduplication (latest version per CASEID), pediatric subsetting, ROR/IC025 |
| `faers_phase2.py`, `faers_phase2_zippass.py`, `faers_phase2_assemble.py` | revision-round sensitivity/context analyses: quarter-resolved disproportionality, masking assessment, age-missingness sensitivity, primary-suspect-only exposure, ADHD comparator class, GI control events |
| `faers_phase3_zippass.py`, `faers_phase3_assemble.py` | notoriety-context comparator panel: stable-warning psychotropic classes (SSRIs, second-generation antipsychotics) vs GI control events, pre/post-2023 fold-change |
| `yrbs_pilot.py`, `yrbs_weightcontrol.py` | YRBS 2025 survey-weighted models |
| `yrbs_extended_covariates.py` | YRBS 2025 extended-covariate sensitivity (bullying, NSSI, binge eating, food insecurity, sexual identity) |
| `extract_cis.py`, `clump_harmonize.py` | GLP1R cis-window extraction (locus ±100 kb, r² < 0.3) and harmonization |
| `run_child_mr.py`, `run_child_mr_developmental.py` | wave-specific childhood cis-MR and fixed-infant-instrument developmental scan |
| `extract_adult_meta_mvp.py`, `run_adult_meta_mvp_mr.py`, `run_adult_relaxed.py` | adult BMI meta-analysis instruments and MR |
| `run_ld_gls.py` | LD-corrected generalized-least-squares sensitivity |
| `run_coloc.py`, `run_coloc_susie.py` | GLP1R ±250-kb colocalization: single-variant coloc.abf posteriors and SuSiE-RSS multi-signal sensitivity |
| `run_steiger_filter.py` | standard Steiger directionality filtering (Fisher z test on variant–trait correlations) |
| `supp_rev3.py` | supplementary table assembly patches (comparator panel, colocalization posteriors) |
| `plot_submission_fig1.py`, `plot_developmental_pub.py`, `plot_submission_figS1.py`, `plot_figure_s2.py`, `plot_graphical_abstract.py` | manuscript and supplement figures plus graphical abstract (vector PDF) |

> **Figure-numbering note:** script filenames predate the final supplement renumbering. In the published supplement, the quarter-resolved FAERS figure (`plot_figure_s2.py`) appears as **Figure S1**, and the wave-specific GLP1R cis-MR forest plot (`plot_submission_figS1.py`) appears as **Figure S2**.

## Data sources

- **FAERS**: FDA public quarterly data files (<https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-public-dashboard>)
- **YRBS 2025**: CDC public-use file (<https://www.cdc.gov/yrbs/>)
- **MoBa childhood BMI GWAS summary statistics**: Helgeland et al., *Nature Metabolism* 2022 (Norwegian Institute of Public Health)
- **iPSYCH suicide-attempt GWAS**: available to qualified researchers on application to the iPSYCH consortium (Danish institutional affiliation and Danish Data Protection Agency approval required)
- **GIANT / Million Veteran Program BMI**: public repositories

Raw GWAS summary statistics and FAERS raw quarterly files are not redistributed; result-level derived quantities are archived on Zenodo (DOI above). SNP-level effect sizes from restricted-access iPSYCH data are excluded per its data-use terms.

## Status

Manuscript under review. Companion dataset: <https://doi.org/10.5281/zenodo.23008872>.
