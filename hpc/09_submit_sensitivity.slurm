#!/bin/bash
#SBATCH --job-name=alemu_sens
#SBATCH --partition=cpu
#SBATCH --array=1-240%120
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/sensitivity_%A_%a.out
#SBATCH --error=logs/sensitivity_%A_%a.err
set -eo pipefail
source /etc/profile
module load R/R-4.6.1
set -u
traits=(DTF DTM HSW PH PPP YPPlnt)
cvs=(CV0 CV1 CV2 CV00)
idx=$((SLURM_ARRAY_TASK_ID-1));partition=$((idx/24+1));within=$((idx%24))
trait=${traits[$((within/4))]};cv=${cvs[$((within%4))]}
: "${PROJECT_DIR:?Export PROJECT_DIR as the reproducibility archive path before sbatch}"
Rscript "$PROJECT_DIR/05_nested_meta_gwas_m4/08_fit_top100_top1000.R" "$partition" "$trait" "$cv"
