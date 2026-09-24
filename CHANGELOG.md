# Changelog

This archive was independently re-executed end-to-end (excluding stages that
require the raw VCF, the raw plot-level phenotype file, or HPC access -- see
`verification/README.md`). Four real bugs were found during that
verification and fixed here. All three prevented the archive from running
its own documented example commands; none change the reported analytical
results, since the manuscript's numbers were generated before these
scripts were archived in their broken state.

## Fixed

### 1. `02_phenotypic_analysis/estimate_variance_components.R`
**Bug**: `getvc()` extracted variance components using
`grepl(pattern, rownames(vc), fixed = TRUE)` with pattern `"Genotype$"`.
With `fixed = TRUE`, this searches for the *literal* substring `"Genotype$"`,
which never appears in ASReml's output (the actual row name is `"Genotype"`,
with no literal dollar sign). This meant `genetic_variance` and
`H2_entry_mean` were always `NA`, regardless of input data.
**Fix**: changed the pattern to `"^Genotype$"` with `fixed = FALSE` (proper
regex, anchored to avoid also matching `Genotype:location` etc).
**Verified**: re-run on both the 209- and 190-genotype datasets now returns
real heritability estimates matching the manuscript's Table 1 values
(209-genotype run) and revealing the correction needed for the 190-genotype
values that should have been reported (see the manuscript's response to
Reviewer 2).

### 2. `03_relationship_matrices/01_build_kernels.R`
**Bug**: crashed with `subscript out of bounds` when run exactly as
documented in this README's own example command. `second_stage_BLUEs_Y.csv`
contains 209 originally-phenotyped accessions; the genome-wide and
haplotype-tagged genotype matrices contain only the 190 that passed
genotype-level QC. The script indexed the genotype matrices directly by
all 209 phenotype-file accessions without first checking they existed in
the genotype panel, so the 19 accessions with phenotypes but no genotype
record (`GG2`, `ICC11879`, `ICC12492`, ... -- the same 19 removed by
`--mind 0.2`) caused an immediate error.
**Fix**: added an explicit intersection step, restricting both the
phenotype data and the accession list to those present in both genotype
panels before any indexing.
**Verified**: now runs to completion, producing an object with the correct
190-accession, 1,710-observation structure, with eigenvalues correlating
>0.998 with the officially supplied `inputs/step5_panels.rds`.

### 3. `05_nested_meta_gwas_m4/01_fit_alemu_m4.R` and `08_fit_top100_top1000.R`
**Bug**: both scripts default `scratch_root` (or `scratch`) to
`file.path("/scratch", Sys.getenv("USER"), "alemu_tmp")` when the
`SLURM_TMPDIR` environment variable is not set. `/scratch/$USER/...` is
specific to the original HPC cluster and does not exist on any other
machine, so `dir.create()` silently fails and the subsequent `setwd()`
throws `cannot change working directory` immediately on the first
environment of the first fold -- before any GWAS computation happens at
all.
**Fix**: fall back to `file.path(tempdir(), "alemu_tmp")` (or
`"alemu_sensitivity_tmp"`), which is portable and works on any machine.
**Verified**: `01_fit_alemu_m4.R` now runs past the previous crash point
into genuine BLINK/GAPIT computation, confirmed via GAPIT's own runtime
diagnostic output and successful completion of at least one full
environment-level GWAS run (see `verification/m4_gwas_test_DTF_CV1_partition1/`).

### 4. `04_cross_validation_and_prediction/01_fit_m1_m3_exact_partitions.R`
**Bug**: `fit_one()` called `BGLR(...)` without a `saveAt` argument. BGLR
defaults to writing its MCMC trace files (`mu.dat`, `varE.dat`,
`ETA_1_varB.dat`, `ETA_2_varB.dat`, ...) to the current working directory
using the same fixed filenames on every call. This script calls `fit_one()`
up to 5 times per fold (M1, M2 x 2 panels, M3 x 2 panels) and is intended to
loop over all 6 traits x 4 CV schemes x 5 folds in one invocation -- every
one of those calls would silently overwrite the previous call's trace files
and litter whatever directory the script happens to be run from. This was
caught because it left `ETA_*.dat`, `mu.dat`, and `varE.dat` files at the
repository root during verification.
**Fix**: `fit_one()` now takes a `tag` argument and gives every BGLR call a
unique `saveAt` prefix under `tempdir()`, deleting the trace files
immediately after each fit (matching the pattern already used correctly in
`05_nested_meta_gwas_m4/01_fit_alemu_m4.R` and `08_fit_top100_top1000.R`).
**Verified**: re-run with the fix produces no stray files at the repository
root and completes fold-by-fold as before.

## Also fixed (test correctness, not a pipeline bug)

### `04_cross_validation_and_prediction/validate_step5.R`
**Bug**: reported `ALL PASS: ...` whenever `output_mechanism/` (a directory
of historical ground-truth output) did not exist locally, because it
silently checked zero files and zero-checked/zero-failed was treated as a
pass.
**Fix**: now raises an error instead of reporting a false pass when no
files were actually checked.

## Not changed

- `01_genotype_qc_and_tagging/` and `02_phenotypic_analysis/01_within_environment_blues.R`
  require external raw data files (the source VCF and the raw plot-level
  phenotype file respectively) that are intentionally not duplicated in
  this archive. See `verification/README.md` for what was independently
  confirmed about these stages despite not being able to execute them
  directly.

## Removed

### `05_nested_meta_gwas_m4/assemble_supplementary_figure_S1.m`
This was not MATLAB -- it was an Objective-C command-line tool using macOS-only
AppKit/PDFKit frameworks to combine the 24 individual two-page Manhattan/Q-Q
PDFs into the single assembled Supplementary Figure S1. It was undocumented
in the "Software" section of this README (which lists only R, BGLR, GAPIT,
ASReml-R, PLINK, Haploview, and Java), would not compile on Linux or Windows,
and nothing else in the repository referenced it. Its only output -- the
assembled Figure S1 PDF -- is already present as a finished deliverable in
`outputs/supplementary_figures/`, so its removal does not reduce what a
reviewer can reproduce or inspect.
