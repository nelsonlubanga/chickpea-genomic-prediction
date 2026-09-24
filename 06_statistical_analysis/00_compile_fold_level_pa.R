#!/usr/bin/env Rscript

## Compile M1-M4 prediction CSVs to one fold-level predictive-ability file.
## Usage:
##   Rscript 00_compile_fold_level_pa.R <M1-M3_root> <M4_project_root> <top_k> <output.csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop(paste(
  "Usage: Rscript 00_compile_fold_level_pa.R",
  "<M1-M3_root> <M4_project_root> <top_k> <output.csv>"))

base_root <- normalizePath(args[[1L]], mustWork = TRUE)
m4_root <- normalizePath(args[[2L]], mustWork = TRUE)
top_k <- as.integer(args[[3L]])
outfile <- args[[4L]]
stopifnot(top_k %in% c(100L, 500L, 1000L))

traits <- c("DTF", "DTM", "HSW", "PH", "PPP", "YPPlnt")
cvs <- c("CV1", "CV2", "CV0", "CV00")
panels <- c(chr = "genome-wide", tags = "haplotype-tagged")

fold_pa <- function(path) {
  if (!file.exists(path)) stop("Missing prediction file: ", path)
  d <- read.csv(path, check.names = FALSE)
  required <- c("environment", "y", "yHat")
  miss <- setdiff(required, names(d))
  if (length(miss)) stop("Missing columns in ", path, ": ", paste(miss, collapse = ", "))
  env_cor <- vapply(split(d, d$environment), function(z) {
    ok <- complete.cases(z$y, z$yHat)
    if (sum(ok) < 5L || sd(z$y[ok]) == 0 || sd(z$yHat[ok]) == 0) return(NA_real_)
    suppressWarnings(cor(z$y[ok], z$yHat[ok]))
  }, numeric(1))
  env_cor <- env_cor[is.finite(env_cor)]
  if (!length(env_cor)) NA_real_ else mean(env_cor)
}

rows <- vector("list", 10L * 6L * 4L * 5L * 8L)
i <- 0L
add <- function(model, panel, trait, cv, repetition, fold, path) {
  i <<- i + 1L
  rows[[i]] <<- data.frame(Model = model, Panel = panel, Trait = trait, CV = cv,
                           Repetition = repetition, Fold = fold, PA = fold_pa(path))
}

for (r in 1:10) for (trait in traits) for (cv in cvs) for (fold in 1:5) {
  bdir <- file.path(base_root, paste0("p", r), trait, cv, paste0("fold_", fold))
  mdir <- file.path(m4_root, "results", sprintf("partition_%02d", r), trait, cv,
                    paste0("fold_", fold))
  for (panel_code in names(panels)) {
    panel <- panels[[panel_code]]
    add("M1", panel, trait, cv, r, fold, file.path(bdir, "M1.csv"))
    add("M2", panel, trait, cv, r, fold,
        file.path(bdir, paste0("M2_", panel_code, ".csv")))
    add("M3", panel, trait, cv, r, fold,
        file.path(bdir, paste0("M3_", panel_code, ".csv")))
    add("M4", panel, trait, cv, r, fold,
        file.path(mdir, sprintf("M4_%s_top%d.csv", panel_code, top_k)))
  }
}

ans <- do.call(rbind, rows[seq_len(i)])
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
write.csv(ans, outfile, row.names = FALSE)
cat("Wrote", nrow(ans), "fold-level rows to", outfile, "\n")
