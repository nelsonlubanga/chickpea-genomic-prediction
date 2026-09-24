#!/bin/bash
#SBATCH --job-name=m4_post
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/postprocess_%j.out
#SBATCH --error=logs/postprocess_%j.err
set -eo pipefail
source /etc/profile
module load R/R-4.6.1
set -u
: "${PROJECT_DIR:?Export PROJECT_DIR as the reproducibility archive path before sbatch}"
Rscript "$PROJECT_DIR/05_nested_meta_gwas_m4/03_gwas_plots_and_reviewer_tables.R"
