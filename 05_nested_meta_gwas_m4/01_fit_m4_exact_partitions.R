#!/usr/bin/env Rscript

## Alemu-inspired M4 extended to multi-environment GxE prediction.
## Usage: Rscript 01_fit_m4_exact_partitions.R <partition 1..10> <trait> <CV>
## Marker discovery is nested within each outer CV training set. Environment-
## specific BLINK results are combined by signed, sample-size-weighted Stouffer
## meta-analysis. The top K markers are fitted as a BRR (shrunk random) term.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript 01_fit_m4_exact_partitions.R <partition> <trait> <CV>")
partition <- as.integer(args[1]); trait <- args[2]; cv_name <- args[3]
traits <- c("DTF", "DTM", "HSW", "PH", "PPP", "YPPlnt")
cv_names <- c("CV0", "CV1", "CV2", "CV00")
if (is.na(partition) || partition < 1L || partition > 10L) stop("partition must be 1..10")
if (!trait %in% traits) stop("Unknown trait: ", trait)
if (!cv_name %in% cv_names) stop("Unknown CV scheme: ", cv_name)

suppressPackageStartupMessages({library(GAPIT); library(BGLR); library(data.table)})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR"), mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
out_root <- file.path(project_dir, "results")
scratch_root <- Sys.getenv("SLURM_TMPDIR", unset = "")
if (!nzchar(scratch_root)) {
  ## Fall back to a portable temp directory when not running under Slurm
  ## (e.g. a reviewer reproducing a single environment/fold locally). The
  ## previous default (/scratch/$USER/alemu_tmp) is HPC-specific and does
  ## not exist on other machines, causing setwd() to fail outright.
  scratch_root <- file.path(tempdir(), "alemu_tmp")
}
dir.create(scratch_root, recursive = TRUE, showWarnings = FALSE)
top_k <- as.integer(Sys.getenv("TOP_K", unset = "500"))
min_proxy_r2 <- as.numeric(Sys.getenv("MIN_PROXY_R2", unset = "0.5"))
ld_window <- as.numeric(Sys.getenv("LD_WINDOW_BP", unset = "500000"))
n_iter <- as.integer(Sys.getenv("N_ITER", unset = "5000"))
burn_in <- as.integer(Sys.getenv("BURN_IN", unset = "2000"))
thin <- as.integer(Sys.getenv("THIN", unset = "5"))
set.seed(12345 + partition)

base_panels <- readRDS(file.path(input_dir, "step5_panels.rds"))
Y <- base_panels$Y; E <- base_panels$E; L <- base_panels$L; panels <- base_panels$panels
accessions <- readRDS(file.path(input_dir, "step5_base.rds"))$accessions
parts12 <- readRDS(file.path(input_dir, "step5_partitions.rds"))[[partition]]
parts00 <- readRDS(file.path(input_dir, "step5_partitions_cv0_cv00.rds"))[[partition]]

GD <- readRDS(file.path(input_dir, "genome_wide_genotypes_190x409128.rds"))
stopifnot(is.data.frame(GD), names(GD)[1] == "Taxa")
GD <- GD[match(accessions, GD$Taxa), , drop = FALSE]
stopifnot(identical(as.character(GD$Taxa), accessions))
X <- as.matrix(GD[-1]); rownames(X) <- accessions; storage.mode(X) <- "double"
for (j in which(colSums(is.na(X)) > 0L)) X[is.na(X[,j]),j] <- mean(X[,j], na.rm=TRUE)
GM_full <- as.data.frame(fread(file.path(input_dir, "genome_wide_marker_map.tsv")))
rownames(GM_full) <- GM_full$SNP
## GAPIT 4.1 requires the numerical-genotype marker map to contain exactly
## marker identifier, chromosome and position. Alleles remain in GM_full for
## reporting but must not be passed as extra GM columns.
GM <- GM_full[,c("SNP","Chromosome","Position")]
tags <- scan(file.path(input_dir, "haplotype_tagged_marker_list.txt"), what="", quiet=TRUE)

fold_mask <- function(fold) {
  if (cv_name == "CV1") return(list(test=which(parts12$CV1 == fold), excluded=integer()))
  if (cv_name == "CV2") return(list(test=which(parts12$CV2 == fold), excluded=integer()))
  if (cv_name == "CV0") return(list(test=which(parts00$CV0 == fold), excluded=integer()))
  list(
    test=which(parts00$geno_fold == fold & parts00$env_fold == fold),
    excluded=which(xor(parts00$geno_fold == fold, parts00$env_fold == fold))
  )
}

run_env_gwas <- function(train_df, env_name, fold) {
  d <- train_df[train_df$env == env_name & !is.na(train_df[[trait]]), c("strain", trait)]
  names(d) <- c("Taxa", trait)
  d <- d[match(intersect(accessions, d$Taxa), d$Taxa), , drop=FALSE]
  if (nrow(d) < 20L) return(NULL)
  gd <- GD[match(d$Taxa, GD$Taxa), , drop=FALSE]
  gdir <- file.path(scratch_root, sprintf("p%02d_%s_%s_f%d", partition, trait, cv_name, fold), env_name)
  dir.create(gdir, recursive=TRUE, showWarnings=FALSE)
  old <- getwd(); on.exit(setwd(old), add=TRUE); setwd(gdir)
  ans <- tryCatch(
    GAPIT(Y=d, GD=gd, GM=GM, PCA.total=3, model="BLINK",
          Random.model=FALSE, Multi_iter=FALSE, file.output=FALSE),
    error=function(e) {message("GAPIT failed for ", env_name, ": ", conditionMessage(e)); NULL})
  if (is.null(ans) || is.null(ans$GWAS) || !nrow(ans$GWAS)) return(NULL)
  z <- as.data.table(ans$GWAS)
  if ("effect" %in% names(z) && !"Effect" %in% names(z)) setnames(z,"effect","Effect")
  if ("maf" %in% names(z) && !"MAF" %in% names(z)) setnames(z,"maf","MAF")
  z$environment <- env_name; z$n <- nrow(d)
  setwd(old)
  unlink(gdir, recursive=TRUE)
  as.data.frame(z)
}

meta_rank <- function(env_results) {
  usable <- Filter(function(z) all(c("SNP", "P.value") %in% names(z)), env_results)
  if (!length(usable)) stop("No usable environment-specific GWAS results")
  pieces <- lapply(usable, function(z) {
    p <- pmax(pmin(as.numeric(z$P.value), 1), .Machine$double.xmin)
    effect <- if ("Effect" %in% names(z)) as.numeric(z$Effect) else rep(1, nrow(z))
    sign_effect <- sign(effect); sign_effect[!is.finite(sign_effect) | sign_effect == 0] <- 1
    data.frame(SNP=z$SNP, z=sign_effect*qnorm(p/2, lower.tail=FALSE), w=sqrt(z$n[1]),
               environment=z$environment[1], P.value=p, Effect=effect)
  })
  long <- rbindlist(pieces)
  meta <- long[, .(Z_meta=sum(w*z, na.rm=TRUE)/sqrt(sum(w^2)),
                         N_environments=.N,
                         min_P=min(P.value, na.rm=TRUE),
                         consistent_sign=length(unique(sign(Effect[is.finite(Effect) & Effect != 0]))) <= 1),
                     by=SNP]
  meta[, P_meta := 2*pnorm(abs(Z_meta), lower.tail=FALSE)]
  setorder(meta, P_meta, -N_environments, min_P, SNP)
  list(ranking=as.data.frame(meta), per_environment=as.data.frame(long))
}

map_to_tags <- function(snps) {
  mapped <- character()
  details <- list()
  for (s in snps) {
    if (s %in% tags) {
      mapped <- c(mapped, s); details[[length(details)+1L]] <- data.frame(discovery=s, proxy=s, r2=1)
      next
    }
    if (!s %in% rownames(GM_full)) next
    chrom <- as.character(GM_full[s,"Chromosome"]); pos <- as.numeric(GM_full[s,"Position"])
    cand <- tags[startsWith(tags, paste0(chrom, "_"))]
    if (!length(cand)) next
    cand_pos <- suppressWarnings(as.numeric(sub("^.*_", "", cand)))
    near <- cand[is.finite(cand_pos) & abs(cand_pos-pos) <= ld_window]
    if (!length(near)) next
    r2 <- vapply(near, function(tg) suppressWarnings(cor(X[,s], X[,tg])^2), numeric(1))
    r2[!is.finite(r2)] <- 0; best <- which.max(r2)
    if (r2[best] >= min_proxy_r2) {
      mapped <- c(mapped, near[best])
      details[[length(details)+1L]] <- data.frame(discovery=s, proxy=near[best], r2=r2[best])
    }
  }
  list(markers=unique(mapped), details=if(length(details)) do.call(rbind, details) else data.frame())
}

fit_m4 <- function(markers, panel_name, y_train, y, test, out_dir) {
  p <- panels[[panel_name]]
  ETA <- list(
    list(X=E, model="BRR"), list(X=L, model="BRR"),
    list(V=p$EVD.G$vectors, d=p$EVD.G$values, model="RKHS"),
    list(V=p$EVD.GxE$vectors, d=p$EVD.GxE$values, model="RKHS")
  )
  if (length(markers)) {
    xs <- X[,markers,drop=FALSE]
    xs <- xs[,apply(xs,2,sd)>0,drop=FALSE]
    if (ncol(xs)) ETA[[length(ETA)+1L]] <- list(X=scale(L %*% xs),model="BRR")
  }
  prefix <- file.path(scratch_root, paste0("bglr_", partition, "_", trait, "_", cv_name, "_", panel_name, "_"))
  fm <- BGLR(y=y_train, ETA=ETA, nIter=n_iter, burnIn=burn_in, thin=thin,
             verbose=FALSE, saveAt=prefix)
  unlink(Sys.glob(paste0(prefix, "*.dat")))
  test <- test[!is.na(y[test])]
  write.csv(data.frame(testing=test, Individual=Y$strain[test], environment=Y$env[test],
                       y=y[test], yHat=fm$yHat[test]),
            file.path(out_dir, paste0("M4_",panel_name,"_top",top_k,".csv")), row.names=FALSE)
}

y <- Y[[trait]]
for (fold in 1:5) {
  out_dir <- file.path(out_root, sprintf("partition_%02d",partition), trait, cv_name, paste0("fold_",fold))
  done <- file.path(out_dir, ".complete")
  if (file.exists(done)) {message("Skipping completed ", out_dir); next}
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  mask <- fold_mask(fold); y_train <- y; y_train[c(mask$test,mask$excluded)] <- NA
  train_rows <- which(!is.na(y_train))
  train_df <- Y[train_rows,c("strain","env",trait)]
  envs <- sort(unique(train_df$env))
  message(sprintf("p%d %s %s fold %d: %d training records, %d environments", partition,trait,cv_name,fold,length(train_rows),length(envs)))
  env_results <- lapply(envs, function(e) run_env_gwas(train_df,e,fold)); names(env_results) <- envs
  combined <- meta_rank(env_results)
  ranking <- combined$ranking
  selected <- head(ranking$SNP, min(top_k,nrow(ranking)))
  write.csv(ranking, file.path(out_dir,"meta_analysis_ranking.csv"), row.names=FALSE)
  write.csv(combined$per_environment[combined$per_environment$SNP %in% selected,],
            file.path(out_dir,"selected_marker_environment_results.csv"),row.names=FALSE)
  writeLines(selected,file.path(out_dir,paste0("selected_top",top_k,".txt")))
  fit_m4(intersect(selected,colnames(X)),"chr",y_train,y,mask$test,out_dir)
  tag_map <- map_to_tags(selected)
  write.csv(tag_map$details,file.path(out_dir,"tag_proxy_mapping.csv"),row.names=FALSE)
  if (length(tag_map$markers)) fit_m4(tag_map$markers,"tags",y_train,y,mask$test,out_dir)
  else {
    fit_m4(character(),"tags",y_train,y,mask$test,out_dir)
    writeLines("No tag proxy met the 500-kb and r2 >= 0.5 criteria; M4 reverted to M3.",
               file.path(out_dir,"M4_tags_M3_fallback.txt"))
  }
  writeLines(c(paste("TOP_K",top_k),"META signed sample-size-weighted Stouffer","SELECTED_EFFECT BRR",
               paste("N_ITER",n_iter),paste("BURN_IN",burn_in),paste("THIN",thin)),
             file.path(out_dir,"model_specification.txt"))
  file.create(done)
}
