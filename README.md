# Reproducible workflow: chickpea multi-environment genomic prediction

Pipeline for chickpea genomic prediction (models M1-M4, GWAS-informed marker
panels) across nine environments, including genotype QC, within-environment
BLUE estimation, and cross-validation.

**This pipeline has been independently re-executed end-to-end** (excluding
stages that require the raw VCF, the raw plot-level phenotype file, or HPC
access). Five bugs that broke or silently corrupted the archive's own
documented example commands were found and fixed in the process. See `CHANGELOG.md` for exactly what
changed and why, and `verification/README.md` for the commands run, the
real output produced, and honest disclosure of what could and could not be
confirmed locally.

## Purpose

This archive contains the final analysis code supporting the manuscript. It
compares a 409,128-SNP genome-wide panel with a 9,845-marker LD-based
haplotype-tagged panel for six traits, four prediction models, and four
cross-validation scenarios in 190 genotypes evaluated in nine environments.

This is the submission workflow for the final **top-500 meta-GWAS/BRR M4
analysis**. It supersedes the earlier FDR/marker-union/fixed-effect M4 workflow.

## Final M4 specification

Within every outer cross-validation training fold:

1. BLINK GWAS is run separately in each available training environment using
   training phenotypes only and the first three marker principal components.
2. Environment-specific evidence is combined by signed, sample-size-weighted
   Stouffer meta-analysis.
3. Markers are ranked by `P_meta`; combined statistics are used for ranking,
   not formal inference.
4. The top 500 markers are fitted as the BRR-shrunk `S` component in
   `E + L + G + GxE + S`.
5. The same rankings are reused for top-100 and top-1,000 sensitivity analyses.
6. For the haplotype-tagged panel, an exact tag is retained; otherwise the
   strongest proxy within 500 kb is used only when `r2 >= 0.50`. Duplicate
   proxies are collapsed.

Validation phenotypes never contribute to GWAS discovery or marker ranking.

## Directory map

| Directory | Contents |
|---|---|
| `01_genotype_qc_and_tagging/` | PLINK QC, two-pass Haploview Tagger selection, and dosage export |
| `02_phenotypic_analysis/` | Within-environment BLUEs and variance-component/heritability analysis |
| `03_relationship_matrices/` | Environment, line, G, and GxE matrices and eigendecompositions |
| `04_cross_validation_and_prediction/` | Exact CV assignments and M1-M3 fitting |
| `05_nested_meta_gwas_m4/` | Final nested BLINK/meta-GWAS M4, sensitivity models, GWAS audit, and completeness checks |
| `06_statistical_analysis/` | Fold-level PA compilation and mixed-model comparisons |
| `07_figures/` | Final manuscript PA figures |
| `08_supplementary_tables/` | Supplementary Tables S1-S5 and S7 builders |
| `09_environment/` | Software versions, R package versions, and requirement checker |
| `hpc/` | Slurm submission templates used on the HPC |
| `inputs/` | Analysis-ready phenotypes, exact folds, kernels, genotype matrices, marker map, and tag list |
| `reference_results/` | Final fold-level PA and summary/statistical outputs used for manuscript reporting |
| `outputs/supplementary_figures/` | Supplementary Figure S1, its caption, and page manifest |
| `verification/` | Independent re-execution checks, captured outputs, and reviewer-requested GWAS diagnostics (see `verification/README.md`) |

See `DATA_DICTIONARY.md` for every abbreviation, code, and core output column.

## Analysis order

### 1. Validate the supplied environment

From this directory:

```bash
Rscript 09_environment/check_requirements.R
Rscript validate_archive.R
```

`check_requirements.R` checks the complete refitting environment and will flag
ASReml-R or GAPIT when they are unavailable. They are not required merely to
regenerate the final plots and supplementary tables from the supplied results.

To regenerate all reporting outputs in one command:

```bash
bash 00_run_reporting_workflow.sh
```

### 2. Rebuild genotype panels from the source VCF (optional)

The source VCF is not duplicated in this archive. To repeat QC and tagging,
provide `Filt.PanCa.vcf.gz`, PLINK 1.90, Java, and Haploview 4.1, then run the
scripts in `01_genotype_qc_and_tagging/`. Final assertions are 190 genotypes,
409,128 genome-wide SNPs, and 9,845 tag SNPs.

### 3. Within-environment BLUE estimation (not executable by this team)

`02_phenotypic_analysis/01_within_environment_blues.R` documents, but was
**not run by this pipeline**, the model used to obtain within-environment
BLUEs:

```bash
Rscript 02_phenotypic_analysis/01_within_environment_blues.R \
  plot_level_phenotypes.csv outputs/within_environment_BLUEs.csv
```

This first-stage model (line fixed, block random, per environment) was run
by **ICRISAT** on their own raw, plot-level field records -- which include
block assignments this team has never had access to -- before the resulting
BLUEs were shared for this study. The script is included so the model
specification is unambiguous and machine-checkable, not because this team
can execute it: `plot_level_phenotypes.csv` (a file with `environment`,
`line`, `block`, and trait columns) does not exist anywhere in this team's
data holdings. The analysis-ready output of that upstream step -- the
second-stage BLUEs actually used for every downstream analysis in this
repository -- is supplied as `inputs/second_stage_BLUEs_Y.csv`.

#### Across-environment variance components and heritability (Table 1)

Fitted to the BLUEs of the 190 genotyped lines with ASReml-R: environment
fixed; line, line x site and line x year random. Requires ASReml-R.

```bash
Rscript 02_phenotypic_analysis/estimate_variance_components.R \
  verification/second_stage_BLUEs_Y_190_verification.csv \
  reference_results/tables/Heritability_190.csv \
  reference_results/tables/Variance_components_REML_190.csv
```

`Variance_components_REML_190.csv` holds the full REML output for each trait
(component, standard error, z-ratio, boundary status, percentage of total
variance, log-likelihood, convergence). `Heritability_190.csv` holds the
entry-mean H² reported in Table 1.

### 4. Rebuild genomic kernels (optional)

```bash
Rscript 03_relationship_matrices/01_build_kernels.R \
  inputs/second_stage_BLUEs_Y.csv \
  inputs/genome_wide_genotypes_190x409128.rds \
  inputs/haplotype_tagged_genotypes_190x9845.rds \
  outputs/step5_panels.rds
```

The exact kernel object used in the reported analysis is supplied as
`inputs/step5_panels.rds`.

### 5. Refit M1-M3

For one replicate partition:

```bash
TRAITS=DTF CV_SCHEMES=CV1 Rscript \
  04_cross_validation_and_prediction/01_fit_m1_m3_exact_partitions.R \
  inputs/step5_panels.rds inputs/step5_partitions.rds \
  inputs/step5_partitions_cv0_cv00.rds 1 outputs/m1_m3
```

Omit `TRAITS` and `CV_SCHEMES` to run all six traits and all four schemes.

### 6. Refit the primary M4 and sensitivity models

The M4 scripts expect `PROJECT_DIR` to contain `inputs/` and `results/`.

```bash
export PROJECT_DIR="$PWD"
export TOP_K=500 MIN_PROXY_R2=0.5 LD_WINDOW_BP=500000
export N_ITER=5000 BURN_IN=2000 THIN=5
Rscript 05_nested_meta_gwas_m4/01_fit_m4_exact_partitions.R 1 DTF CV1
Rscript 05_nested_meta_gwas_m4/08_fit_top100_top1000.R 1 DTF CV1
```

The complete analysis is 10 repetitions x 6 traits x 4 CV schemes; each task
fits five folds. Adapt the templates in `hpc/` to the local cluster project
path and R module name.

### 7. Regenerate the final figures

```bash
Rscript 07_figures/01_make_final_pa_plots.R
```

This uses the included primary top-500 fold-level file and writes PNG and PDF
versions of the four CV-specific plots and the combined plot to
`reference_results/figures/`. Optional input and output paths can be supplied
as the first and second command-line arguments.

Supplementary Figure S1 contains representative fold-level meta-GWAS
Manhattan and Q-Q plots for repetition 1, fold 1 across all six traits and four
CV schemes. The fixed selection rule, caption, and page order are documented in
`outputs/supplementary_figures/README.md`. The 24 individual per-trait/CV
source PDFs used to assemble it are not included in this repository (each is
~19 MB; ~456 MB total) -- only the compact, already-assembled
`Supplementary_Figure_S1_representative_GWAS_plots.pdf` is kept here. The
complete set of 1,200 fold-level plot PDFs (all partitions, all folds;
approximately 22 GB) is not deposited. It can be regenerated from the
fold-level meta-analysis rankings with
`05_nested_meta_gwas_m4/03_gwas_plots_and_reviewer_tables.R` after running
step 6, and is available from the authors on request.

### 8. Regenerate supplementary tables

The submission workbooks are built from the final fold-level PA, statistical
outputs, and reviewer-facing meta-GWAS files. See
`08_supplementary_tables/01_build_supplementary_tables_S1_S5.R` and the README
inside the deposited supplementary-table directory.

Table S6 (variance components) is produced by
`02_phenotypic_analysis/estimate_variance_components.R` (see step 3). Table S7
(representation of the M4 markers in the haplotype-tagged panel) is rebuilt
from the deposited top-500 selections with:

```bash
Rscript 08_supplementary_tables/02_build_table_S7_marker_composition.R
```

which reproduces the submitted counts exactly (~20 minutes).

### 9. Fold-level meta-GWAS results and Supplementary Figure S2

`05_nested_meta_gwas_m4/10_regenerate_fold_gwas_plots.R` reruns the GWAS
stage of M4 for one fold (no BGLR fits) and writes the full marker ranking,
its Manhattan/Q-Q plots and a check against the deposited top 500. All 1,200
folds were regenerated this way; the rankings and plots are deposited at
Zenodo (see the Licence and citation section), and 1,192 of 1,200 folds
reproduce the original top-500 marker set exactly.
`11_redraw_fold_gwas_plots.R` redraws a fold's plots from its saved ranking.

Supplementary Figure S2 (environment-specific Q-Q plots before aggregation,
and the effect of between-environment correlation on the combined statistic)
is produced by `verification/environment_qq_and_correlation_check.R`.

To compile new prediction files into the fold-level PA structure used by all
reporting scripts, run:

```bash
Rscript 06_statistical_analysis/00_compile_fold_level_pa.R \
  outputs/m1_m3 . 500 outputs/PA_fold_level_M1_M4_top500.csv
```

## Fixed parameters and seeds

- Final genotype QC: chromosomes Ca1-Ca8; `--mind 0.2 --geno 0.1 --maf 0.05`.
- LD tagging: two-pass Haploview Tagger; pairwise `r2 >= 0.50`, with 1-kb
  preliminary thinning.
- BLINK: GAPIT 4.1.0; first three marker PCs; minimum 20 observations per
  training-environment GWAS.
- Meta-analysis: signed, sample-size-weighted Stouffer statistic.
- Primary selected set: top 500; sensitivity sets: top 100 and top 1,000.
- Tag proxy: strongest proxy within 500 kb and `r2 >= 0.50`.
- BGLR: 5,000 iterations, 2,000 burn-in, thinning interval 5.
- M4 seed: `12345 + repetition`; sensitivity seed: `22345 + repetition`.
- Exact reported cross-validation assignments are supplied in `inputs/`.
- PA: within-environment Pearson correlation where at least five validation
  observations are available, averaged equally across eligible environments.
- Statistical model: `PA ~ Model + (1 | Repetition/Fold)` or
  `PA ~ Panel + (1 | Repetition/Fold)`; Tukey adjustment for model contrasts.

## Software

Core software versions are recorded in `09_environment/`. The final workflow
used R 4.6.1, BGLR 1.1.4, GAPIT 4.1.0, ASReml-R 4.2.0.482, PLINK 1.90b6.21,
Haploview 4.1, and Java 19.0.1. ASReml-R is proprietary.

## Reproducibility boundaries

- The raw VCF is not duplicated here. The whole-genome sequence data were
  obtained from Varshney et al. (2021) and are available from that study.
- The raw, plot-level phenotype file (with block assignments) is not simply
  omitted from this archive -- this team never held it. ICRISAT ran the
  within-environment BLUE model (line fixed, block random; see
  `02_phenotypic_analysis/01_within_environment_blues.R` and the manuscript
  Methods) on their own raw plot-level records and shared only the resulting
  BLUEs, supplied here as `inputs/second_stage_BLUEs_Y.csv`. That first-stage
  step is therefore documented in this repository but cannot be independently
  re-executed by anyone working from this archive alone.
- Exact CV assignment objects are supplied because the complete historical
  random-number stream used to generate all ten partitions was not preserved.
- The environment-level selected-marker evidence has 4,860,000 rows and is
  deposited as TSV rather than Excel because it exceeds Excel's row limit.
- Large stochastic refits may show negligible Monte Carlo variation, while the
  included fold assignments and reference results preserve the reported design
  and numerical outputs.
- Two files in `reference_results/` are large (`Supplementary_Table_S5_GWAS_results_underlying_M4.xlsx`,
  ~64 MB, and `gwas/selected_markers_all_folds.tsv`, ~75 MB); both are kept
  here (under GitHub's 100 MB hard limit) and are also included in the
  Figshare snapshot of this repository cited below.
- See `verification/README.md` for exactly which parts of this pipeline have
  been independently re-executed and confirmed, and which parts (genotype QC
  from the raw VCF, within-environment BLUEs from the raw plot-level file,
  and the full HPC-scale CV/GWAS grid) could not be tested locally and why.

## Output-to-manuscript mapping

| Manuscript item | Source |
|---|---|
| Figures 2-5 | `07_figures/01_make_final_pa_plots.R` and top-500 fold-level PA |
| Supplementary Figure S1 | Representative repetition-1/fold-1 fold-level meta-GWAS Manhattan and Q-Q plots |
| Table S1 | Top-500 fold-level and summary PA |
| Table S2 | Top-100 sensitivity PA |
| Table S3 | Top-1,000 sensitivity PA |
| Table S4 | Primary top-500 mixed-model and panel comparisons |
| Table S5 | Fold-specific top-500 meta-GWAS rankings and selection frequencies |

## Licence and citation

Please cite the associated manuscript and the Figshare snapshot of this
repository: `https://doi.org/10.6084/m9.figshare.33311715`.

The complete fold-level meta-GWAS results (full marker ranking and
Manhattan/Q-Q plots for all 1,200 folds) are deposited at Zenodo:
[Zenodo DOI].
