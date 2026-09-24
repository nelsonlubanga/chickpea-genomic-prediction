suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 01_build_kernels.R <phenotypes.csv> <genome_wide.rds> <tagged.rds> <output.rds>")
}
y <- fread(args[1])
if (!all(c("strain", "env") %in% names(y))) stop("Phenotype file needs strain and env columns")

read_panel <- function(path) {
  d <- readRDS(path)
  if (is.data.frame(d) && names(d)[1] == "Taxa") {
    ids <- as.character(d[[1]]); x <- as.matrix(d[-1])
  } else {
    x <- as.matrix(d); ids <- rownames(x)
  }
  rownames(x) <- ids
  storage.mode(x) <- "double"
  for (j in seq_len(ncol(x))) {
    mu <- mean(x[, j], na.rm = TRUE)
    x[is.na(x[, j]), j] <- mu
    s <- sd(x[, j])
    x[, j] <- if (!is.finite(s) || s == 0) 0 else (x[, j] - mu) / s
  }
  x[, apply(x, 2, sd) > 0, drop = FALSE]
}

X_chr_full <- read_panel(args[2])
X_tags_full <- read_panel(args[3])

## The phenotype file (second_stage_BLUEs_Y.csv) contains all originally
## phenotyped accessions; the genotype panels contain only those that passed
## genotype-level QC. Restrict to the intersection so downstream indexing
## does not fail on accessions with phenotypes but no genotype record.
genotyped_accessions <- intersect(rownames(X_chr_full), rownames(X_tags_full))
y <- y[as.character(y$strain) %in% genotyped_accessions]
accessions <- unique(as.character(y$strain))
if (!length(accessions)) stop("No accessions in the phenotype file overlap the supplied genotype panels")

X_chr <- X_chr_full[accessions, , drop = FALSE]
X_tags <- X_tags_full[accessions, , drop = FALSE]

L <- model.matrix(~ factor(strain, levels = accessions) - 1, data = y)
E <- model.matrix(~ factor(env) - 1, data = y)
ZE <- tcrossprod(E)

make_panel <- function(x) {
  G_geno <- tcrossprod(x) / ncol(x)
  G <- tcrossprod(tcrossprod(L, G_geno), L)
  GxE <- G * ZE
  list(G_genotype = G_geno, EVD.G = eigen(G, symmetric = TRUE),
       EVD.GxE = eigen(GxE, symmetric = TRUE))
}

ans <- list(Y = as.data.frame(y), E = E, L = L,
            panels = list(chr = make_panel(X_chr), tags = make_panel(X_tags)))
saveRDS(ans, args[4])
cat("Saved E, L, G and GxE objects for", nrow(y), "observations and", length(accessions), "accessions.\n")
