#!/usr/bin/env Rscript

## Supplementary Table S2: top-500 M4 markers classified as direct, proxy
## (best tag within 500 kb, r2 >= 0.5) or excluded, as in map_to_tags() of M4.
## Usage: Rscript 08_supplementary_tables/02_build_table_S2_marker_composition.R [out.xlsx]

suppressPackageStartupMessages({library(data.table); library(writexl)})
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args)) args[1] else
  "reference_results/supplementary_tables/Supplementary_Table_S2_M4_tagged_panel_marker_composition.xlsx"
ld_window <- 500000; min_proxy_r2 <- 0.5

sel <- fread("reference_results/gwas/selected_markers_all_folds.tsv")
tags <- scan("inputs/haplotype_tagged_marker_list.txt", what = "", quiet = TRUE)
GM <- fread("inputs/genome_wide_marker_map.tsv")
accessions <- readRDS("inputs/step5_base.rds")$accessions
GD <- readRDS("inputs/genome_wide_genotypes_190x409128.rds")
GD <- GD[match(accessions, GD$Taxa), , drop = FALSE]

snps <- unique(sel$SNP)
need <- unique(c(snps, tags))
X <- as.matrix(GD[, need[need %in% names(GD)], drop = FALSE]); storage.mode(X) <- "double"
for (j in which(colSums(is.na(X)) > 0L)) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
rm(GD); invisible(gc())

tag_chr <- sub("_.*$", "", tags)
tag_pos <- suppressWarnings(as.numeric(sub("^.*_", "", tags)))
map <- GM[match(snps, SNP), .(SNP, Chromosome, Position)]

classify <- function(s, chrom, pos) {
  if (s %in% tags) return(list("direct", s, 1))
  if (is.na(chrom)) return(list("excluded", NA_character_, NA_real_))
  near <- tags[tag_chr == chrom & is.finite(tag_pos) & abs(tag_pos - pos) <= ld_window]
  if (!length(near)) return(list("excluded", NA_character_, NA_real_))
  r2 <- vapply(near, function(tg) suppressWarnings(cor(X[, s], X[, tg])^2), numeric(1))
  r2[!is.finite(r2)] <- 0; b <- which.max(r2)
  if (r2[b] >= min_proxy_r2) list("proxy", near[b], unname(r2[b])) else list("excluded", near[b], unname(r2[b]))
}
cls <- map[, {r <- classify(SNP, Chromosome, Position); .(status = r[[1]], proxy = r[[2]], r2 = r[[3]])}, by = SNP]

sel <- merge(sel[, .(Partition, Trait, CV, Fold, SNP)], cls[, .(SNP, status)], by = "SNP")
cv_order <- c("CV0", "CV00", "CV1", "CV2")
by_cv <- dcast(sel[, .N, by = .(Trait, CV, status)], Trait + CV ~ status, value.var = "N", fill = 0L)
by_cv[, total := direct + excluded + proxy]
by_cv[, `:=`(direct_per_fold = round(direct / 50, 1), proxy_per_fold = round(proxy / 50, 1),
             excluded_per_fold = round(excluded / 50, 1))]
by_part <- dcast(sel[, .N, by = .(Trait, CV, Partition, status)], Trait + CV + Partition ~ status,
                 value.var = "N", fill = 0L)
by_part[, total := direct + proxy + excluded]
setcolorder(by_part, c("Trait", "CV", "Partition", "direct", "proxy", "excluded", "total"))
setorder(by_cv, Trait, CV); setorder(by_part, Trait, CV, Partition)
unique_cls <- cls[, .(SNP, status, proxy, r2 = round(r2, 4))]

readme <- data.table(
  Item = c("Title", "Description", "Classification rule", "Columns (By_trait_CV)", "Unique_SNP_classification", "Code"),
  Details = c(
    "Supplementary Table S2. Composition of the GWAS-ranked marker component (S) of model M4 in the haplotype-tagged panel, by trait and cross-validation scheme",
    sprintf("Each of the %s fold-level marker selections of the primary top-500 analysis (10 repetitions x 6 traits x 4 CV schemes x 5 folds x 500 markers) is classified against the %s-marker haplotype-tagged panel.",
            format(nrow(sel), big.mark = ","), format(length(tags), big.mark = ",")),
    "direct: the selected SNP is itself a tag SNP. proxy: otherwise, the strongest tag SNP on the same chromosome within 500 kb has r2 >= 0.5 with it (genotype dosages of the 190 genotypes). excluded: no tag SNP qualifies. Identical to map_to_tags() in 05_nested_meta_gwas_m4/01_fit_m4_exact_partitions.R.",
    "direct, proxy, excluded and total are counts over 10 repetitions x 5 folds (50 folds of 500 markers); *_per_fold are the means per fold.",
    "One row per unique selected SNP: status, the best tag SNP within 500 kb (proxy) and its r2.",
    "08_supplementary_tables/02_build_table_S2_marker_composition.R in https://github.com/nelsonlubanga/chickpea-genomic-prediction"))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_xlsx(list(README = readme, By_trait_CV = by_cv, By_trait_CV_partition = by_part,
                Unique_SNP_classification = unique_cls), out_file)
cat("Wrote", out_file, "\n"); print(by_cv)
