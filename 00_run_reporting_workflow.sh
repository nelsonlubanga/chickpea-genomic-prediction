#!/usr/bin/env bash
set -euo pipefail

archive_dir="$(cd "$(dirname "$0")" && pwd)"
cd "$archive_dir"

Rscript validate_archive.R
Rscript 06_statistical_analysis/01_test_models_and_panels.R 500 \
  reference_results/data/PA_fold_level_M1_M4_top500.csv \
  outputs/Table_M4_top500_model_panel_comparison.xlsx
Rscript 07_figures/01_make_final_pa_plots.R \
  reference_results/data/PA_fold_level_M1_M4_top500.csv \
  outputs/figures
Rscript 08_supplementary_tables/01_build_supplementary_tables_S1_S5.R \
  "$archive_dir" "$archive_dir/outputs/supplementary_tables" \
  "$archive_dir/reference_results/gwas"

echo "Reporting workflow complete: $archive_dir/outputs"
