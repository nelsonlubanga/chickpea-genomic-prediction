#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
results_dir <- file.path(project, "results")
map_file <- file.path(project, "inputs", "genome_wide_marker_map.tsv")
out_dir <- file.path(project, "reviewer_outputs", "Supplementary_Figure_S1_sources_with_thresholds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("DTF", "DTM", "HSW", "PH", "PPP", "YPPlnt")
cvs <- c("CV1", "CV2", "CV0", "CV00")
partition_to_show <- 1L
fold_to_show <- 1L
chr_order <- paste0("Ca", 1:8)
chr_cols <- rep(c("#2C7FB8", "#7FCDBB"), 4)

map <- fread(map_file)
ranking_files <- list.files(
  results_dir,
  pattern = "meta_analysis_ranking\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

parse_key <- function(f) {
  x <- strsplit(sub(paste0("^", results_dir, "/?"), "", f), "/", fixed = TRUE)[[1]]
  list(
    Partition = as.integer(sub("partition_", "", x[1])),
    Trait = x[2],
    CV = x[3],
    Fold = as.integer(sub("fold_", "", x[4]))
  )
}

keys <- lapply(ranking_files, parse_key)
keep <- vapply(keys, function(k) {
  k$Partition == partition_to_show && k$Fold == fold_to_show &&
    k$Trait %in% traits && k$CV %in% cvs
}, logical(1))
ranking_files <- ranking_files[keep]
keys <- keys[keep]

if (length(ranking_files) != length(traits) * length(cvs)) {
  stop("Expected 24 ranking files for partition 1/fold 1; found ", length(ranking_files))
}

manifest <- vector("list", length(ranking_files))
for (i in seq_along(ranking_files)) {
  f <- ranking_files[i]
  key <- keys[[i]]
  d <- merge(fread(f), map, by = "SNP", all.x = TRUE, sort = FALSE)
  pdat <- d[Chromosome %in% chr_order & is.finite(Position) & is.finite(P_meta)]
  pdat[, Chromosome := factor(Chromosome, levels = chr_order)]
  chr_len <- pdat[, .(len = max(Position, na.rm = TRUE)), by = Chromosome]
  chr_len[, offset := shift(cumsum(len), fill = 0)]
  pdat <- merge(pdat, chr_len[, .(Chromosome, offset)], by = "Chromosome", sort = FALSE)
  pdat[, x := Position + offset]

  display_trait <- if (key$Trait == "YPPlnt") "YPP" else key$Trait
  out_file <- file.path(
    out_dir,
    sprintf("%s_%s_partition01_fold1.pdf", key$Trait, key$CV)
  )
  pdf(out_file, width = 12, height = 6, onefile = TRUE)
  plot(
    pdat$x, -log10(pmax(pdat$P_meta, .Machine$double.xmin)),
    pch = 20, cex = 0.28, col = chr_cols[as.integer(pdat$Chromosome)],
    xaxt = "n", xlab = "Chromosome", ylab = expression(-log[10](italic(P))),
    main = sprintf("%s %s partition 1 fold 1: meta-GWAS", display_trait, key$CV)
  )
  axis(1, at = pdat[, mean(range(x)), by = Chromosome]$V1, labels = chr_order)
  d[, FDR_meta := p.adjust(P_meta, method = "BH")]
  bonf_p <- 0.05 / sum(is.finite(d$P_meta))
  abline(h = -log10(bonf_p), col = "#D7301F", lty = 2, lwd = 1.5)
  bh_pass <- d[is.finite(P_meta) & FDR_meta <= 0.05, P_meta]
  if (length(bh_pass)) {
    abline(h = -log10(max(bh_pass)), col = "#762A83", lty = 3, lwd = 1.5)
  }
  legend(
    "topright",
    legend = c("Bonferroni 0.05", if (length(bh_pass)) "BH-FDR 5%"),
    col = c("#D7301F", if (length(bh_pass)) "#762A83"),
    lty = c(2, if (length(bh_pass)) 3), lwd = 1.5, bty = "n", cex = 0.8
  )
  obs <- sort(pmax(pmin(d$P_meta, 1), .Machine$double.xmin))
  exp <- ppoints(length(obs))
  plot(
    -log10(exp), -log10(obs), pch = 20, cex = 0.35,
    xlab = expression(Expected~~-log[10](italic(P))),
    ylab = expression(Observed~~-log[10](italic(P))),
    main = sprintf("%s %s partition 1 fold 1: Q-Q", display_trait, key$CV)
  )
  abline(0, 1, col = "red", lwd = 1.5)
  dev.off()

  manifest[[i]] <- data.table(
    Trait = display_trait, CV = key$CV, Repetition = 1L, Fold = 1L,
    Source = f, Figure_page = basename(out_file)
  )
}

manifest <- rbindlist(manifest)
manifest[, Trait := factor(Trait, levels = c("DTF", "DTM", "HSW", "PH", "PPP", "YPP"))]
manifest[, CV := factor(CV, levels = c("CV1", "CV2", "CV0", "CV00"))]
setorder(manifest, Trait, CV)
manifest[, `:=`(Trait = as.character(Trait), CV = as.character(CV))]
fwrite(manifest, file.path(out_dir, "Supplementary_Figure_S1_manifest.tsv"), sep = "\t")
cat("Created", nrow(manifest), "representative Manhattan/Q-Q pages in", out_dir, "\n")
