suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 02_export_imputed_matrices.R <work_dir> <output_dir>")
}
work_dir <- normalizePath(args[1], mustWork = TRUE)
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_plink_raw <- function(path) {
  x <- fread(path)
  taxa <- x$IID
  x[, c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE") := NULL]
  setnames(x, sub("_[ACGT]$", "", names(x)))
  for (j in seq_len(ncol(x))) {
    set(x, j = j, value = as.numeric(x[[j]]))
    if (anyNA(x[[j]])) set(x, which(is.na(x[[j]])), j, mean(x[[j]], na.rm = TRUE))
  }
  data.frame(Taxa = taxa, x, check.names = FALSE)
}

gw <- read_plink_raw(file.path(work_dir, "ca1_8_qc_mind02_numeric.raw"))
tags <- read_plink_raw(file.path(work_dir, "final_tagged_ca1_8_r05_v2_numeric.raw"))
stopifnot(nrow(gw) == 190L, ncol(gw) - 1L == 409128L, !anyNA(gw))
stopifnot(nrow(tags) == 190L, ncol(tags) - 1L == 9845L, !anyNA(tags))
stopifnot(identical(as.character(gw$Taxa), as.character(tags$Taxa)))

saveRDS(gw, file.path(out_dir, "genome_wide_genotypes_190x409128.rds"))
saveRDS(tags, file.path(out_dir, "haplotype_tagged_genotypes_190x9845.rds"))
cat("Wrote mean-imputed matrices to", out_dir, "\n")
