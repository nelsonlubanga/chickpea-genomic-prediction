#!/usr/bin/env Rscript

## Reviewer 2, Point 5: environment-specific (pre-aggregation) Q-Q plots, and
## the effect of correlation among environments on the combined statistic.
##
## Input: per-environment BLINK results for one trait on the full data
## (written by correlation_aware_reranking_check.R; same GAPIT call as M4).
## For each SNP the M4 statistic is Z = sum(w z) / sqrt(sum(w^2)), w = sqrt(n).
## Allowing for a correlation matrix R among environment-level z-scores gives
## Z_R = sum(w z) / sqrt(w' R w). Every SNP is tested in the same environments
## with the same weights, so Z_R = c * Z for a single constant c: the ranking,
## and hence any top-k marker set, cannot change; only the calibration of P does.
## R is estimated from markers with |z| < 2 in every environment (approximately
## null markers); the phenotypic correlation between environments is shown for
## comparison.
##
## Usage: Rscript environment_qq_and_correlation_check.R [trait] [project_dir]
## Outputs (verification/environment_qq_<trait>/):
##   Supplementary_Figure_S2_environment_QQ_<trait>.pdf / .png
##   environment_qq_summary.csv, combined_statistic_summary.csv

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
trait <- if (length(args) >= 1L) args[[1L]] else "DTF"
project_dir <- if (length(args) >= 2L) args[[2L]] else normalizePath(".", mustWork = TRUE)
in_dir <- file.path(project_dir, "verification", sprintf("correlation_aware_reranking_%s", trait))
out_dir <- file.path(project_dir, "verification", sprintf("environment_qq_%s", trait))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

env_results <- readRDS(file.path(in_dir, "raw_env_gwas_results.rds"))
R_pheno <- as.matrix(fread(file.path(in_dir, "environment_phenotypic_correlation.csv")), rownames = 1)

## Environment-level signed z-scores, as in meta_rank() of the M4 script.
envs <- names(env_results)
snps <- sort(env_results[[1]]$SNP)
Z <- sapply(env_results, function(z) {
  z <- as.data.table(z)[match(snps, SNP)]
  p <- pmax(pmin(as.numeric(z$P.value), 1), .Machine$double.xmin)
  s <- sign(as.numeric(z$Effect)); s[!is.finite(s) | s == 0] <- 1
  s * qnorm(p / 2, lower.tail = FALSE)
})
stopifnot(!anyNA(Z))
n <- sapply(env_results, function(z) z$n[1])
w <- sqrt(n)
lambda_gc <- function(z) median(z^2) / qchisq(0.5, 1)

null <- apply(abs(Z), 1, max) < 2
R_null <- cor(Z[null, ])
combine <- function(R) {
  zc <- as.vector(Z %*% w) / sqrt(drop(t(w) %*% R %*% w))
  p <- 2 * pnorm(-abs(zc))
  list(z = zc, p = p)
}
naive <- combine(diag(length(envs)))
adj_null <- combine(R_null)
adj_pheno <- combine(R_pheno[envs, envs])
off <- function(R) mean(R[upper.tri(R)])

env_summary <- data.table(trait, environment = envs, n = n,
                          lambda_GC = round(apply(Z, 2, lambda_gc), 3))
comb_summary <- data.table(
  trait,
  statistic = c("naive Stouffer (used in M4)", "correlation-adjusted: null-marker z correlation",
                "correlation-adjusted: phenotypic correlation"),
  mean_between_env_correlation = round(c(0, off(R_null), off(R_pheno[envs, envs])), 3),
  scale_factor_vs_naive = round(c(1, sqrt(sum(w^2) / drop(t(w) %*% R_null %*% w)),
                                  sqrt(sum(w^2) / drop(t(w) %*% R_pheno[envs, envs] %*% w))), 4),
  lambda_GC = round(sapply(list(naive, adj_null, adj_pheno), function(x) lambda_gc(x$z)), 3),
  n_P_below_5e8 = sapply(list(naive, adj_null, adj_pheno), function(x) sum(x$p < 5e-8)),
  spearman_rank_vs_naive = sapply(list(naive, adj_null, adj_pheno),
                                  function(x) cor(rank(-abs(naive$z)), rank(-abs(x$z)), method = "spearman")),
  top500_overlap_vs_naive = sapply(list(naive, adj_null, adj_pheno), function(x)
    length(intersect(order(naive$p)[1:500], order(x$p)[1:500])))
)
fwrite(env_summary, file.path(out_dir, "environment_qq_summary.csv"))
fwrite(comb_summary, file.path(out_dir, "combined_statistic_summary.csv"))
print(env_summary); print(comb_summary)

## Figure: nine environment-specific Q-Q plots, then the combined statistic
## before and after the correlation adjustment.
qq_panel <- function(p, main, sub) {
  obs <- -log10(sort(pmax(pmin(p, 1), .Machine$double.xmin)))
  exp <- -log10(ppoints(length(obs)))
  plot(exp, obs, pch = 20, cex = 0.3, xlab = expression(Expected~~-log[10](italic(P))),
       ylab = expression(Observed~~-log[10](italic(P))))
  title(main, line = 1.6, cex.main = 1)
  abline(0, 1, col = "red", lwd = 1.5)
  mtext(sub, side = 3, line = 0.3, cex = 0.75)
}
draw <- function() {
  par(mfrow = c(3, 4), mar = c(4.2, 4.2, 3.4, 1))
  for (e in envs) {
    qq_panel(2 * pnorm(-abs(Z[, e])), sprintf("%s: %s", trait, e),
             bquote(list(italic(n) == .(unname(n[e])), lambda[GC] == .(sprintf("%.2f", lambda_gc(Z[, e]))))))
  }
  qq_panel(naive$p, "Combined (naive Stouffer, as in M4)", bquote(lambda[GC] == .(sprintf("%.2f", lambda_gc(naive$z)))))
  qq_panel(adj_null$p, "Combined (correlation-adjusted)",
           bquote(lambda[GC] == .(sprintf("%.2f", lambda_gc(adj_null$z))) * ";  marker ranking unchanged"))
  plot.new()
  text(0, 0.5, adj = 0, cex = 0.85, paste(
    "Panels 1-9: single-environment BLINK GWAS",
    "(before aggregation), all genotyped",
    "accessions with a record in the environment.",
    "Panel 10: combined statistic used to rank",
    "markers in M4. Panel 11: the same statistic",
    "allowing for the correlation among",
    sprintf("environment-level z-scores (mean r = %.2f).", off(R_null)), sep = "\n"))
}
stem <- file.path(out_dir, sprintf("Supplementary_Figure_S2_environment_QQ_%s", trait))
png(paste0(stem, ".png"), width = 16, height = 12, units = "in", res = 200); draw(); dev.off()
pdf(paste0(stem, ".pdf"), width = 16, height = 12)
par(mar = c(0, 0, 0, 0)); plot.new()
rasterImage(png::readPNG(paste0(stem, ".png")), 0, 0, 1, 1, interpolate = FALSE)
dev.off()
cat("Done. Output in", out_dir, "\n")
