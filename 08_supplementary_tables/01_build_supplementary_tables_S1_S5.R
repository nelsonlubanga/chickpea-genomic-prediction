#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1L) args[[1L]] else ".", mustWork = TRUE)
out_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(root, "outputs", "supplementary_tables")
src_dir <- if (length(args) >= 3L) args[[3L]] else
  file.path(root, "reference_results", "gwas")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readme_frame <- function(title, description, rows) {
  data.frame(
    Item = c("Title", "Description", names(rows)),
    Details = c(title, description, unname(rows)),
    check.names = FALSE
  )
}

pretty_data <- function(d) {
  d <- as.data.frame(d)
  if ("Trait" %in% names(d)) d$Trait[d$Trait == "YPPlnt"] <- "YPP"
  if ("Panel" %in% names(d)) {
    d$Panel[d$Panel == "genome-wide"] <- "Genome-wide"
    d$Panel[d$Panel == "haplotype-tagged"] <- "Haplotype-tagged"
  }
  d
}

write_pa_table <- function(k, number) {
  means_file <- file.path(root, "reference_results", "tables",
                          sprintf("Means_M4_top%d_by_CV_trait_panel_model.csv", k))
  folds_file <- file.path(root, "reference_results", "data",
                          sprintf("PA_fold_level_M1_M4_top%d.csv", k))
  means <- pretty_data(fread(means_file))
  folds <- pretty_data(fread(folds_file))
  names(means)[names(means) == "PA"] <- "Mean_predictive_ability"
  names(means)[names(means) == "SD"] <- "SD_predictive_ability"
  names(means)[names(means) == "N"] <- "N_folds"
  names(folds)[names(folds) == "PA"] <- "Predictive_ability"
  readme <- readme_frame(
    sprintf("Supplementary Table S%d. Predictive abilities from models M1-M4 using the top %d GWAS-ranked markers in M4", number, k),
    sprintf("Trait-, cross-validation-, panel-, and model-specific predictive abilities for the M4-%d analysis.", k),
    c(
      `Marker specification` = sprintf("M4 contains the %d highest-ranked fold-specific meta-GWAS markers; M1-M3 are unchanged.", k),
      `Marker discovery` = "BLINK GWAS was nested within every training fold and run separately by training environment using training phenotypes only and three marker principal components.",
      `Across-environment ranking` = "Signed, sample-size-weighted Stouffer meta-analysis; combined statistics were used for ranking rather than formal inference.",
      `Predictive ability` = "Within-environment Pearson correlation between observed and predicted phenotypes for environments with at least five validation observations, averaged equally within each fold.",
      `Summary sheet` = "Arithmetic mean, standard deviation, and number of fold-level predictive abilities across 10 repetitions x 5 folds.",
      `Fold_level sheet` = "Complete fold-level predictive abilities used to produce the summary and statistical comparisons.",
      `CV0` = "Tested genotypes in untested environments.",
      `CV1` = "Untested genotypes in observed environments.",
      `CV2` = "Sparse testing: tested genotypes with missing records in observed environments.",
      `CV00` = "Untested genotypes in untested environments.",
      `Models` = "M1: baseline; M2: genomic main effect; M3: genomic main effect plus GxE; M4: M3 plus a BRR-shrunk component for the selected GWAS-ranked markers.",
      `Trait label` = "YPP corresponds to the internal analysis label YPPlnt."
    )
  )
  outfile <- file.path(out_dir, sprintf("Supplementary_Table_S%d_top%d_markers.xlsx", number, k))
  write_xlsx(list(README = readme, Summary = means, Fold_level = folds), outfile, format_headers = TRUE)
  outfile
}

created <- c(
  write_pa_table(500, 1),
  write_pa_table(100, 2),
  write_pa_table(1000, 3)
)

## Preserve the already formatted and statistically validated top-500 workbook.
s4_source <- file.path(root, "reference_results", "tables",
                       "Table_M4_top500_model_panel_comparison.xlsx")
s4_target <- file.path(out_dir, "Supplementary_Table_S4_statistical_tests_top500_markers.xlsx")
if (!file.copy(s4_source, s4_target, overwrite = TRUE)) stop("Could not create Table S4")
created <- c(created, s4_target)

selected_file <- file.path(src_dir, "selected_markers_all_folds.tsv")
frequency_file <- file.path(src_dir, "marker_selection_frequency_across_partitions.tsv")
audit_file <- file.path(src_dir, "GWAS_reproducibility_audit.tsv")
stopifnot(file.exists(selected_file), file.exists(frequency_file), file.exists(audit_file))

selected <- fread(selected_file)
frequency <- fread(frequency_file)
audit <- fread(audit_file)
selected$Trait[selected$Trait == "YPPlnt"] <- "YPP"
frequency$Trait[frequency$Trait == "YPPlnt"] <- "YPP"

s5_readme <- readme_frame(
  "Supplementary Table S5. Genome-wide association results underlying the GWAS-assisted model (M4)",
  "Fold-nested meta-GWAS results for the top 500 markers selected in each trait, cross-validation scheme, repetition, and fold, together with marker-selection frequencies and the reproducibility specification.",
  c(
    `Meta_GWAS_selected sheet` = "The 500 highest-ranked genome-wide SNPs from each of 1,200 analyses (10 repetitions x 6 traits x 4 CV schemes x 5 folds); 600,000 rows in total.",
    `Selection_frequency sheet` = "Frequency with which each SNP was selected across the 50 repetition-fold combinations for each trait and CV scheme.",
    `Reproducibility_audit sheet` = "Software and analytical specifications used for GWAS discovery, meta-analysis, leakage prevention, marker selection, shrinkage, and tag-proxy mapping.",
    `Leakage prevention` = "Every GWAS used training-fold phenotypes only; validation phenotypes were excluded before marker discovery and ranking.",
    `Z_meta` = "Signed, sample-size-weighted Stouffer meta-analysis statistic.",
    `P_meta` = "Two-sided probability calculated from Z_meta and used to rank markers; it was not treated as formal inference because correlations among environment-specific statistics were not explicitly modelled.",
    `N_environments` = "Number of training environments contributing association evidence for the marker.",
    `min_P` = "Smallest environment-specific GWAS P-value for the marker.",
    `consistent_sign` = "TRUE when all available non-zero environment-specific marker-effect estimates had the same sign.",
    `FDR_meta` = "Benjamini-Hochberg adjusted P_meta within the fold-specific ranking; supplied descriptively and not used to select the fixed top-500 set.",
    `Rank` = "Position in the fold-specific meta-GWAS ranking.",
    `Selection_frequency` = "Selections divided by 50 repetition-fold combinations for the relevant trait and CV scheme.",
    `Environment-level evidence` = "The complete 4,860,000-row environment-level file remains archived on the HPC as selected_marker_environment_evidence.tsv because it exceeds the Excel worksheet row limit."
  )
)

s5_target <- file.path(out_dir, "Supplementary_Table_S5_GWAS_results_underlying_M4.xlsx")
write_xlsx(
  list(
    README = s5_readme,
    Meta_GWAS_selected = as.data.frame(selected),
    Selection_frequency = as.data.frame(frequency),
    Reproducibility_audit = as.data.frame(audit)
  ),
  s5_target,
  format_headers = TRUE
)
created <- c(created, s5_target)

cat(paste(created, collapse = "\n"), "\n")
