#!/usr/bin/env Rscript

required <- c(
  "README.md", "DATA_DICTIONARY.md",
  "inputs/second_stage_BLUEs_Y.csv", "inputs/step5_base.rds",
  "inputs/step5_panels.rds", "inputs/step5_partitions.rds",
  "inputs/step5_partitions_cv0_cv00.rds",
  "inputs/genome_wide_genotypes_190x409128.rds",
  "inputs/haplotype_tagged_genotypes_190x9845.rds",
  "inputs/genome_wide_marker_map.tsv",
  "inputs/haplotype_tagged_marker_list.txt",
  "reference_results/data/PA_fold_level_M1_M4_top100.csv",
  "reference_results/data/PA_fold_level_M1_M4_top500.csv",
  "reference_results/data/PA_fold_level_M1_M4_top1000.csv",
  "07_figures/01_make_final_pa_plots.R"
)

missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing required files:\n", paste(missing, collapse = "\n"))

pa <- read.csv("reference_results/data/PA_fold_level_M1_M4_top500.csv")
stopifnot(nrow(pa) == 9600L)
stopifnot(identical(sort(unique(pa$Model)), c("M1", "M2", "M3", "M4")))
stopifnot(identical(sort(unique(pa$CV)), c("CV0", "CV00", "CV1", "CV2")))
stopifnot(length(unique(pa$Repetition)) == 10L, length(unique(pa$Fold)) == 5L)

tags <- scan("inputs/haplotype_tagged_marker_list.txt", what = "", quiet = TRUE)
map_n <- length(readLines("inputs/genome_wide_marker_map.tsv")) - 1L
stopifnot(length(tags) == 9845L, map_n == 409128L)

cat("Archive validation passed.\n")
cat("Top-500 PA rows:", nrow(pa), "\n")
cat("Genome-wide markers:", map_n, "\n")
cat("Haplotype-tagged markers:", length(tags), "\n")
