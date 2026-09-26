#!/usr/bin/env Rscript

## Produces a genuine, pre-aggregation (single-environment) GWAS QQ plot,
## using the EXACT SAME GAPIT/BLINK call as the production pipeline
## (05_nested_meta_gwas_m4/01_fit_m4_exact_partitions.R's run_env_gwas()):
## PCA.total=3, model="BLINK", Random.model=FALSE. No CV-fold masking is
## applied here -- this uses all 190 genotyped accessions with a phenotype
## record in the chosen environment, to demonstrate what a single
## environment's association evidence looks like BEFORE the Stouffer
## meta-analysis combination step that produces the combined P_meta values
## shown in Supplementary Figure S1.
##
## Usage: Rscript run_single_environment_gwas_for_qq.R <trait> <environment>

suppressPackageStartupMessages({library(GAPIT); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
trait <- if (length(args) >= 1L) args[[1L]] else "DTF"
env_name <- if (length(args) >= 2L) args[[2L]] else "2015.Amlaha"
project_dir <- if (length(args) >= 3L) args[[3L]] else normalizePath(".", mustWork = TRUE)
stopifnot(dir.exists(file.path(project_dir, "inputs")))
input_dir <- file.path(project_dir, "inputs")
out_dir <- file.path(project_dir, "verification", sprintf("single_env_gwas_%s_%s", trait, env_name))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project dir:", project_dir, "\n")
cat("Trait:", trait, " Environment:", env_name, "\n")

Y <- fread(file.path(input_dir, "second_stage_BLUEs_Y.csv"))
GD_full <- readRDS(file.path(input_dir, "genome_wide_genotypes_190x409128.rds"))
GM_full <- fread(file.path(input_dir, "genome_wide_marker_map.tsv"))
GM <- GM_full[, .(SNP, Chromosome, Position)]
accessions <- as.character(GD_full$Taxa)

d <- Y[env == env_name & is.finite(get(trait)), .(Taxa = strain, val = get(trait))]
d <- d[match(intersect(accessions, d$Taxa), d$Taxa)]
setnames(d, "val", trait)
cat("n accessions with phenotype in this environment:", nrow(d), "\n")
stopifnot(nrow(d) >= 20L)

gd <- GD_full[match(d$Taxa, GD_full$Taxa), , drop = FALSE]
stopifnot(identical(as.character(gd$Taxa), as.character(d$Taxa)))

old_wd <- getwd()
gapit_dir <- file.path(tempdir(), sprintf("gapit_single_env_%s_%s", trait, env_name))
dir.create(gapit_dir, recursive = TRUE, showWarnings = FALSE)
setwd(gapit_dir)
cat("Running GAPIT/BLINK (this is the slow step, ~1.5-2 min)...\n")
t0 <- Sys.time()
ans <- GAPIT(Y = as.data.frame(d), GD = gd, GM = as.data.frame(GM), PCA.total = 3,
             model = "BLINK", Random.model = FALSE, Multi_iter = FALSE, file.output = FALSE)
cat("GAPIT runtime:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec\n")
setwd(old_wd)

gwas <- as.data.table(ans$GWAS)
if ("effect" %in% names(gwas) && !"Effect" %in% names(gwas)) setnames(gwas, "effect", "Effect")
fwrite(gwas, file.path(out_dir, "single_environment_gwas_full_results.csv"))
cat("Saved full single-environment GWAS results (", nrow(gwas), "markers ) to\n  ",
    file.path(out_dir, "single_environment_gwas_full_results.csv"), "\n")

## QQ plot
obs <- sort(pmax(pmin(gwas$P.value, 1), .Machine$double.xmin))
exp <- ppoints(length(obs))
png(file.path(out_dir, "single_environment_qq_plot.png"), width = 900, height = 900, res = 150)
plot(-log10(exp), -log10(obs), pch = 20, cex = 0.4,
     xlab = expression(Expected~~-log[10](italic(P))),
     ylab = expression(Observed~~-log[10](italic(P))),
     main = sprintf("%s, %s: single-environment Q-Q\n(pre-aggregation, n=%d)", trait, env_name, nrow(d)))
abline(0, 1, col = "red", lwd = 1.5)
dev.off()
cat("Saved QQ plot to", file.path(out_dir, "single_environment_qq_plot.png"), "\n")

## Genomic inflation factor, for a concrete, quotable number
chisq <- qchisq(1 - obs, df = 1)
lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
cat("Genomic inflation factor (lambda_GC):", round(lambda_gc, 3), "\n")
writeLines(sprintf("trait=%s environment=%s n=%d n_markers=%d lambda_GC=%.3f",
                    trait, env_name, nrow(d), nrow(gwas), lambda_gc),
           file.path(out_dir, "summary.txt"))
cat("Done.\n")
