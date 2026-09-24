#!/bin/bash
#SBATCH --job-name=alemu_check
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/completeness_%j.out
#SBATCH --error=logs/completeness_%j.err
set -eo pipefail
source /etc/profile
module load R/R-4.6.1
set -u
: "${PROJECT_DIR:?Export PROJECT_DIR as the reproducibility archive path before sbatch}"
Rscript "$PROJECT_DIR/05_nested_meta_gwas_m4/05_refresh_completeness.R"
