#!/usr/bin/env bash
# Reproduce the final 409,128-SNP and 9,845-tag panels used in the revision.
# Required environment variables:
#   RAW_VCF, PLINK, JAVA, HAPLOVIEW_JAR
# Optional: WORK_DIR (default: ./work), KEEP_FILE (default: keep209.txt)
set -euo pipefail

: "${RAW_VCF:?Set RAW_VCF to Filt.PanCa.vcf.gz}"
: "${PLINK:?Set PLINK to PLINK v1.90b6.21}"
: "${JAVA:=java}"
: "${HAPLOVIEW_JAR:?Set HAPLOVIEW_JAR to Haploview.jar v4.1}"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WORK_DIR=${WORK_DIR:-"${SCRIPT_DIR}/work"}
KEEP_FILE=${KEEP_FILE:-"${SCRIPT_DIR}/keep209.txt"}
CHROMS=(CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8)
R2=0.5
BP_SPACE=1000
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

"$PLINK" --allow-extra-chr --const-fid --vcf "$RAW_VCF" --vcf-half-call m \
  --keep "$KEEP_FILE" --make-bed --out panel209_raw

"$PLINK" --allow-extra-chr --bfile panel209_raw \
  --chr CA1 CA2 CA3 CA4 CA5 CA6 CA7 CA8 \
  --mind 0.2 --geno 0.1 --maf 0.05 \
  --make-bed --out ca1_8_qc_mind02

test "$(wc -l < ca1_8_qc_mind02.fam | tr -d ' ')" = "190"
test "$(wc -l < ca1_8_qc_mind02.bim | tr -d ' ')" = "409128"

"$PLINK" --allow-extra-chr --bfile ca1_8_qc_mind02 \
  --recode A --out ca1_8_qc_mind02_numeric

mkdir -p thin_1kb tagger_pass1
for chr in "${CHROMS[@]}"; do
  "$PLINK" --allow-extra-chr --bfile ca1_8_qc_mind02 --chr "$chr" \
    --bp-space "$BP_SPACE" --output-missing-phenotype 0 \
    --recode --out "thin_1kb/${chr}"
  awk '{print $2, $4}' "thin_1kb/${chr}.map" > "thin_1kb/${chr}.info"
  "$JAVA" -jar "$HAPLOVIEW_JAR" -n \
    -pedfile "thin_1kb/${chr}.ped" -info "thin_1kb/${chr}.info" \
    -pairwiseTagging -tagrsqcutoff "$R2" -out "tagger_pass1/${chr}"
done
cat tagger_pass1/*.TESTS | sort -u > candidate_tags_r05.txt

"$PLINK" --allow-extra-chr --bfile ca1_8_qc_mind02 \
  --extract candidate_tags_r05.txt --make-bed --out candidate_panel_r05

mkdir -p candidate_by_chr tagger_pass2
for chr in "${CHROMS[@]}"; do
  "$PLINK" --allow-extra-chr --bfile candidate_panel_r05 --chr "$chr" \
    --output-missing-phenotype 0 --recode --out "candidate_by_chr/${chr}"
  awk '{print $2, $4}' "candidate_by_chr/${chr}.map" > "candidate_by_chr/${chr}.info"
  "$JAVA" -jar "$HAPLOVIEW_JAR" -n \
    -pedfile "candidate_by_chr/${chr}.ped" -info "candidate_by_chr/${chr}.info" \
    -pairwiseTagging -tagrsqcutoff "$R2" -out "tagger_pass2/${chr}"
done
cat tagger_pass2/*.TESTS | sort -u > final_tags_r05_v2.txt

"$PLINK" --allow-extra-chr --bfile ca1_8_qc_mind02 \
  --extract final_tags_r05_v2.txt --make-bed --out final_tagged_ca1_8_r05_v2
"$PLINK" --allow-extra-chr --bfile final_tagged_ca1_8_r05_v2 \
  --recode A --out final_tagged_ca1_8_r05_v2_numeric

test "$(wc -l < final_tagged_ca1_8_r05_v2.bim | tr -d ' ')" = "9845"
echo "Verified final panels: 190 accessions; 409,128 genome-wide SNPs; 9,845 tag SNPs."
