# Independent re-execution verification

This directory documents an independent, from-scratch re-execution of every
stage of this pipeline that does not require the raw VCF, the raw plot-level
phenotype file, or an HPC cluster. It exists so a reader can confirm the
code actually does what the manuscript and README claim, without having to
run the full 10-repetition x 6-trait x 4-CV-scheme x 5-fold grid themselves.

Five real bugs were found and fixed in the course of this verification
(see `../CHANGELOG.md` for details). Every result below reflects the
**fixed** code, re-run and captured fresh.

## What was verified, and how to reproduce each check

### 1. Environment and archive integrity

```bash
Rscript 09_environment/check_requirements.R   # all 9 packages detected, exact versions
Rscript validate_archive.R                    # 9,600 PA rows; 409,128 / 9,845 marker counts confirmed
```
Both pass cleanly with no modification.

### 2. Genomic relationship kernels (`03_relationship_matrices/01_build_kernels.R`)

The script as originally archived crashed on its own documented example
command (`subscript out of bounds`), because it indexed the 190-accession
genotype matrix using all 209 originally-phenotyped accessions without first
restricting to the intersection. Fixed (see CHANGELOG). Re-run:

```bash
Rscript 03_relationship_matrices/01_build_kernels.R \
  inputs/second_stage_BLUEs_Y.csv \
  inputs/genome_wide_genotypes_190x409128.rds \
  inputs/haplotype_tagged_genotypes_190x9845.rds \
  <output.rds>
```

Result: completes and produces an object with the correct dimensions
(1,710 observations, 190 accessions), matching `inputs/step5_panels.rds`
exactly in structure. The rebuilt eigenvalues correlate **0.9999**
(genome-wide panel) and **0.9982** (haplotype-tagged panel) with the
official `step5_panels.rds`, tracking the same ranking but not bit-identical
(official values run ~3-15% larger) -- most likely a difference in the exact
marker-standardization convention between this archived script and whatever
originally produced `step5_panels.rds`. This is disclosed here rather than
hidden: the two are extremely close but not a byte-perfect reproduction.

### 3. M1-M3 cross-validated prediction (`04_cross_validation_and_prediction/01_fit_m1_m3_exact_partitions.R`)

Full re-run for one trait/CV scheme/partition:

```bash
TRAITS=DTF CV_SCHEMES=CV1 Rscript 04_cross_validation_and_prediction/01_fit_m1_m3_exact_partitions.R \
  inputs/step5_panels.rds inputs/step5_partitions.rds inputs/step5_partitions_cv0_cv00.rds \
  1 <output_dir>
```

Result: completed all 5 folds x 5 model/panel files (25 files total; see
`m1_m3_test_DTF_CV1_partition1/`). Predictive ability (Pearson correlation,
M3, genome-wide panel, DTF, CV1) recomputed directly from the raw
observed/predicted columns in these files:

| Fold | n  | PA (M3, genome-wide) |
|------|----|----------------------|
| 1    | 341 | 0.856 |
| 2    | 340 | 0.808 |
| 3    | 342 | 0.786 |
| 4    | 339 | 0.821 |
| 5    | 341 | 0.819 |

These values fall within the plausible range described in the manuscript
for a highly heritable trait under CV1 (untested genotypes, observed
environments).

A second bug was found while inspecting this test's output: `fit_one()`
called `BGLR()` without a `saveAt` prefix, so its MCMC trace files
(`mu.dat`, `varE.dat`, `ETA_*.dat`) were written to the current working
directory using fixed names on every call -- 5 calls per fold, colliding
and overwriting each other, and littering the repository root (this is how
it was caught). Fixed by giving every call a unique `saveAt` prefix under
`tempdir()` (see CHANGELOG). Re-verified after the fix: no stray files at
the repository root, and the script completes fold-by-fold exactly as
before.

### 4. Nested meta-GWAS / M4 (`05_nested_meta_gwas_m4/01_fit_m4_exact_partitions.R`)

The script as originally archived hardcoded `/scratch/$USER/m4_tmp` as
its temp-directory fallback, which does not exist off the original HPC
cluster and causes an immediate `cannot change working directory` error.
Fixed to fall back to `tempdir()` (see CHANGELOG). Re-run:

```bash
export PROJECT_DIR=<project_dir_with_inputs_and_results>
export TOP_K=500 MIN_PROXY_R2=0.5 LD_WINDOW_BP=500000
export N_ITER=500 BURN_IN=200 THIN=5   # reduced from the production 5000/2000/5 to keep the spot-check tractable
Rscript 05_nested_meta_gwas_m4/01_fit_m4_exact_partitions.R 1 DTF CV1
```

Result: progresses correctly through real BLINK/GAPIT computation --
confirmed the runtime's own diagnostic line ("There are 190 common
individuals in genotype, phenotype and CV files") on every environment
processed, and completed at least one full environment's GAPIT run
end-to-end ("GAPIT has done all analysis!!!") before this spot-check was
stopped. See `m4_gwas_test_DTF_CV1_partition1/gapit_blink_progress.log`
for the captured console output.

**Runtime note**: a single trait/CV/partition combination took
approximately 1.5-2 minutes per training environment (45 environment-level
GWAS runs for the full 5-fold loop, i.e. roughly 60-75 minutes locally on a
single machine). The complete study is 10 partitions x 6 traits x 4 CV
schemes = 240 such combinations, which is why the production analysis
requires HPC array parallelization (see `hpc/`) rather than sequential
local execution -- this is not a flaw, just a scale reality worth having a
concrete number for.

### 5. Reporting stage (figures, statistics, supplementary tables)

```bash
bash 00_run_reporting_workflow.sh
```

All four steps (`validate_archive.R`, `06_statistical_analysis/01_test_models_and_panels.R`,
`07_figures/01_make_final_pa_plots.R`, `08_supplementary_tables/01_build_supplementary_tables_S1_S5.R`)
run to completion using only the files already in `reference_results/` and
`inputs/` -- no external data or HPC access needed. Output file sizes match
the deposited versions.

### 6. Single-environment (pre-aggregation) GWAS Q-Q check

```bash
Rscript verification/run_single_environment_gwas_for_qq.R DTF 2015.Amlaha
```

Runs the same GAPIT/BLINK call as production (`PCA.total = 3`,
`model = "BLINK"`, `Random.model = FALSE`) for one environment, without
CV-fold masking, to show what one environment's association evidence looks
like before the Stouffer combination used in Supplementary Figure S1.

Result (`single_env_gwas_DTF_2015.Amlaha/`): n = 188 accessions, 409,128
markers, genomic inflation factor **lambda_GC = 1.086**. Full results and the
Q-Q plot are saved alongside `summary.txt`.

### 7. Correlation-aware vs naive Stouffer re-ranking

```bash
Rscript verification/correlation_aware_reranking_check.R DTF .
```

Runs BLINK in each of the nine DTF environments on full (non-CV-masked)
data. It then compares the production naive Stouffer statistic
(denominator `sqrt(sum(w^2))`) with a correlation-aware version
(denominator `sqrt(w' R w)`), where `R` is the between-environment
phenotypic correlation matrix (`environment_phenotypic_correlation.csv`).

Result (`correlation_aware_reranking_DTF/`): top-500 overlap 500/500,
Spearman rank correlation 1.0000.

**Limitation:** this agreement holds by construction, so it is not
evidence of robustness. Every one of the 409,128 markers is tested in all
nine environments with the same per-environment weights. The
correlation-aware denominator is therefore the same constant for every
marker, and `Z_corr` is exactly `0.4742 x Z_naive` throughout
(checked directly from `naive_vs_correlation_aware_ranking.csv`). This
uniform rescaling deflates the combined Z, so it changes the absolute
`P_meta` values, but it cannot change the ranking. A check that could
change the ranking would need marker-specific correlation structure, for
example correlations of the test statistics estimated from null markers,
or a multivariate or random-effects meta-analysis.

## What was not independently verified, and why

- **`01_genotype_qc_and_tagging/`**: requires the raw VCF (`Filt.PanCa.vcf.gz`),
  which is intentionally not duplicated in this archive (see main README,
  "Reproducibility boundaries"). The PLINK log from the original QC run
  (`--chr CA1-CA8 --mind 0.2 --geno 0.1 --maf 0.05`, 19 people removed,
  409,128 markers and 190 people remaining) is consistent with every number
  reported in the manuscript and this archive.
- **`02_phenotypic_analysis/01_within_environment_blues.R`**: this script
  documents, but was never executed by, this pipeline. **ICRISAT** ran this
  exact model (`fixed = y ~ line, random = ~ block`, per environment) on
  their own raw, plot-level field records -- which include block assignments
  -- and shared only the resulting BLUEs (`inputs/second_stage_BLUEs_Y.csv`)
  with this study. This was confirmed directly: an exhaustive search (18,794
  CSV and 56 XLSX files across the full project tree) found zero files with
  a `block` column anywhere, and the second-stage BLUEs were shown to be
  byte-identical, to full floating-point precision, to values in
  `3000_Phenotyping_data_QC_locations updated_3_31_2025.xlsx` -- a
  genotype-level (already-aggregated, no plot/block structure) file, which
  is itself consistent with ICRISAT having already applied the block
  adjustment before this file was produced. The model specification is
  internally consistent with the manuscript's Methods description (once
  correctly attributed to ICRISAT), but the underlying raw plot-level
  records cannot be re-verified from anything in this team's possession.
- **The complete 10x6x4x5 M1-M4/GWAS grid**: confirmed tractable per-cell
  (see above) but requires HPC-scale parallel compute to complete in full,
  exactly as the main README states.
