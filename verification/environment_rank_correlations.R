#!/usr/bin/env Rscript

## Spearman rank correlations of genotype BLUEs between pairs of environments,
## by trait (Discussion: G×E for DTF and HSW).
## Usage: Rscript verification/environment_rank_correlations.R

Y <- readRDS("inputs/step5_panels.rds")$Y
traits <- c("DTF", "DTM", "HSW", "PH", "PPP", "YPPlnt")
res <- do.call(rbind, lapply(traits, function(tr) {
  w <- tapply(Y[[tr]], list(Y$strain, Y$env), mean)
  r <- cor(w, use = "pairwise.complete.obs", method = "spearman")
  v <- r[upper.tri(r)]
  data.frame(trait = tr, pairs = length(v), mean = round(mean(v), 2),
             min = round(min(v), 2), max = round(max(v), 2), negative = sum(v < 0))
}))
print(res, row.names = FALSE)
