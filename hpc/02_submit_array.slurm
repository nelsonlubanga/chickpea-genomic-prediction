#!/bin/bash
#SBATCH --job-name=alemu_m4
#SBATCH --partition=cpu
#SBATCH --array=1-240%24
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/m4_%A_%a.out
#SBATCH --error=logs/m4_%A_%a.err

set -eo pipefail
source /etc/profile
module load R/R-4.6.1
set -u

traits=(DTF DTM HSW PH PPP YPPlnt)
cvs=(CV0 CV1 CV2 CV00)
idx=$((SLURM_ARRAY_TASK_ID-1))
partition=$((idx/24+1))
within=$((idx%24))
trait=${traits[$((within/4))]}
cv=${cvs[$((within%4))]}

: "${PROJECT_DIR:?Export PROJECT_DIR as the reproducibility archive path before sbatch}"
export TOP_K=500
export MIN_PROXY_R2=0.5
export LD_WINDOW_BP=500000
export N_ITER=5000
export BURN_IN=2000
export THIN=5

Rscript "$PROJECT_DIR/05_nested_meta_gwas_m4/01_fit_m4_exact_partitions.R" "$partition" "$trait" "$cv"
