#!/usr/bin/env Rscript

## Reviewer Point 5, third ask: does accounting for correlations among
## environment-specific association statistics change marker ranking,
## relative to the production pipeline's naive (independence-assuming)
## signed sample-size-weighted Stouffer meta-analysis?
##
## Method: run the exact same per-environment BLINK GWAS as production
## (PCA.total=3, model="BLINK", Random.model=FALSE) on the full (non-CV-
## masked) training data for every environment with a record for one
## representative trait. Compute the naive Stouffer Z_meta (as in
## 05_nested_meta_gwas_m4/01_fit_m4_exact_partitions.R's meta_rank()), and a
## correlation-aware version that uses the phenotypic correlation between
## environments (estimated from the shared genotypes' real BLUEs) as a
## proxy for the correlation between their association test statistics --
## a standard approach when independent-sample GWAS replicates are not
## available to estimate that correlation directly. Compare the resulting
## top-500 marker rankings.
##
## Usage: Rscript correlation_aware_reranking_check.R <trait> <project_dir>

suppressPackageStartupMessages({library(GAPIT); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
trait <- if (length(args) >= 1L) args[[1L]] else "DTF"
project_dir <- if (length(args) >= 2L) args[[2L]] else normalizePath(".", mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
out_dir <- file.path(project_dir, "verification", sprintf("correlation_aware_reranking_%s", trait))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

Y <- fread(file.path(input_dir, "second_stage_BLUEs_Y.csv"))
GD_full <- readRDS(file.path(input_dir, "genome_wide_genotypes_190x409128.rds"))
GM_full <- fread(file.path(input_dir, "genome_wide_marker_map.tsv"))
GM <- GM_full[, .(SNP, Chromosome, Position)]
accessions <- as.character(GD_full$Taxa)

envs <- sort(unique(Y[is.finite(get(trait)), env]))
cat("Trait:", trait, "| environments:", paste(envs, collapse=", "), "\n")

## --- Step 1: phenotypic correlation between environments (real data) ---
wide <- dcast(Y[strain %in% accessions], strain ~ env, value.var = trait)
env_mat <- as.matrix(wide[, ..envs])
rownames(env_mat) <- wide$strain
phen_cor <- cor(env_mat, use = "pairwise.complete.obs")
cat("\nPhenotypic correlation between environments (real data, pairwise complete):\n")
print(round(phen_cor, 2))
fwrite(as.data.table(phen_cor, keep.rownames = "env"), file.path(out_dir, "environment_phenotypic_correlation.csv"))

## --- Step 2: run per-environment GWAS on full (non-CV-masked) data ---
run_one <- function(env_name) {
  d <- Y[env == env_name & is.finite(get(trait)), .(Taxa = strain, val = get(trait))]
  d <- d[match(intersect(accessions, d$Taxa), d$Taxa)]
  setnames(d, "val", trait)
  if (nrow(d) < 20L) return(NULL)
  gd <- GD_full[match(d$Taxa, GD_full$Taxa), , drop = FALSE]
  gdir <- file.path(tempdir(), sprintf("gapit_full_%s_%s", trait, env_name))
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(gdir)
  cat(sprintf("[%s] Running GAPIT/BLINK, n=%d...\n", env_name, nrow(d)))
  ans <- tryCatch(
    GAPIT(Y = as.data.frame(d), GD = gd, GM = as.data.frame(GM), PCA.total = 3,
          model = "BLINK", Random.model = FALSE, Multi_iter = FALSE, file.output = FALSE),
    error = function(e) { message("GAPIT failed for ", env_name, ": ", conditionMessage(e)); NULL })
  setwd(old); unlink(gdir, recursive = TRUE)
  if (is.null(ans) || is.null(ans$GWAS) || !nrow(ans$GWAS)) return(NULL)
  z <- as.data.table(ans$GWAS)
  if ("effect" %in% names(z) && !"Effect" %in% names(z)) setnames(z, "effect", "Effect")
  z$environment <- env_name; z$n <- nrow(d)
  z
}

env_results <- lapply(envs, run_one)
names(env_results) <- envs
env_results <- Filter(Negate(is.null), env_results)
cat("\nCompleted GWAS for", length(env_results), "of", length(envs), "environments\n")
saveRDS(env_results, file.path(out_dir, "raw_env_gwas_results.rds"))

## --- Step 3: naive Stouffer (production method, verbatim from meta_rank()) ---
pieces <- lapply(env_results, function(z) {
  p <- pmax(pmin(as.numeric(z$P.value), 1), .Machine$double.xmin)
  effect <- if ("Effect" %in% names(z)) as.numeric(z$Effect) else rep(1, nrow(z))
  sign_effect <- sign(effect); sign_effect[!is.finite(sign_effect) | sign_effect == 0] <- 1
  data.table(SNP = z$SNP, environment = z$environment[1],
             z = sign_effect * qnorm(p / 2, lower.tail = FALSE), w = sqrt(z$n[1]))
})
long <- rbindlist(pieces)

naive <- long[, .(Z_naive = sum(w * z, na.rm = TRUE) / sqrt(sum(w^2)), N_env = .N), by = SNP]
naive[, P_naive := 2 * pnorm(abs(Z_naive), lower.tail = FALSE)]
setorder(naive, P_naive)
naive[, rank_naive := .I]

## --- Step 4: correlation-aware combination ---
## For a SNP observed in environments e1..ek (weights w, z-scores z), the
## naive method divides by sqrt(sum(w^2)) (assumes Cov(z_i,z_j)=0, i!=j).
## The correlation-aware version divides by sqrt(w' R w), where R is the
## environment correlation matrix estimated in Step 1, restricted to the
## environments that actually contributed to this SNP.
present_envs <- long[, .(envs = list(environment)), by = SNP]
corr_z <- function(snp_envs, snp_z, snp_w) {
  R <- phen_cor[snp_envs, snp_envs, drop = FALSE]
  denom <- sqrt(as.numeric(t(snp_w) %*% R %*% snp_w))
  sum(snp_w * snp_z) / denom
}
long_by_snp <- split(long, long$SNP)
Z_corr <- vapply(long_by_snp, function(d) corr_z(d$environment, d$z, d$w), numeric(1))
corr_dt <- data.table(SNP = names(Z_corr), Z_corr = Z_corr)
corr_dt[, P_corr := 2 * pnorm(abs(Z_corr), lower.tail = FALSE)]
setorder(corr_dt, P_corr)
corr_dt[, rank_corr := .I]

## --- Step 5: compare rankings ---
merged <- merge(naive[, .(SNP, P_naive, rank_naive)], corr_dt[, .(SNP, P_corr, rank_corr)], by = "SNP")
top500_naive <- merged[rank_naive <= 500, SNP]
top500_corr <- merged[rank_corr <= 500, SNP]
overlap <- length(intersect(top500_naive, top500_corr))
rank_cor_spearman <- cor(merged$rank_naive, merged$rank_corr, method = "spearman")

fwrite(merged, file.path(out_dir, "naive_vs_correlation_aware_ranking.csv"))
summary_txt <- c(
  sprintf("Trait: %s", trait),
  sprintf("Environments used: %s", paste(names(env_results), collapse = ", ")),
  sprintf("Total markers ranked: %d", nrow(merged)),
  sprintf("Top-500 overlap between naive and correlation-aware ranking: %d / 500 (%.1f%%)",
          overlap, 100 * overlap / 500),
  sprintf("Spearman rank correlation (naive vs correlation-aware), all markers: %.4f", rank_cor_spearman)
)
writeLines(summary_txt, file.path(out_dir, "summary.txt"))
cat("\n", paste(summary_txt, collapse = "\n"), "\n", sep = "")
cat("\nDone. Output in", out_dir, "\n")
