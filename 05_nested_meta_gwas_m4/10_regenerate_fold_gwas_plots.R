#!/usr/bin/env Rscript

## Rerun the GWAS stage of M4 (no model fits) for one fold: ranking, plots and
## a check of the top 500 against selected_markers_all_folds.tsv.
## Usage: Rscript 10_regenerate_fold_gwas_plots.R <partition> <trait> <CV> <fold>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript 10_regenerate_fold_gwas_plots.R <partition> <trait> <CV> <fold>")
partition <- as.integer(args[1]); trait <- args[2]; cv_name <- args[3]; fold <- as.integer(args[4])
stopifnot(partition %in% 1:10, trait %in% c("DTF","DTM","HSW","PH","PPP","YPPlnt"),
          cv_name %in% c("CV0","CV1","CV2","CV00"), fold %in% 1:5)

suppressPackageStartupMessages({library(GAPIT); library(data.table)})

project_dir <- normalizePath(Sys.getenv("PROJECT_DIR", unset = "."), mustWork = TRUE)
input_dir <- file.path(project_dir, "inputs")
out_dir <- file.path(project_dir, "regenerated_gwas",
                     sprintf("partition_%02d", partition), trait, cv_name, paste0("fold_", fold))
done <- file.path(out_dir, ".complete")
if (file.exists(done)) quit(save = "no")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
scratch_root <- file.path(tempdir(), "gwas_tmp")
dir.create(scratch_root, recursive = TRUE, showWarnings = FALSE)

base_panels <- readRDS(file.path(input_dir, "step5_panels.rds"))
Y <- base_panels$Y
accessions <- readRDS(file.path(input_dir, "step5_base.rds"))$accessions
parts12 <- readRDS(file.path(input_dir, "step5_partitions.rds"))[[partition]]
parts00 <- readRDS(file.path(input_dir, "step5_partitions_cv0_cv00.rds"))[[partition]]
GD <- readRDS(file.path(input_dir, "genome_wide_genotypes_190x409128.rds"))
GD <- GD[match(accessions, GD$Taxa), , drop = FALSE]
stopifnot(identical(as.character(GD$Taxa), accessions))
GM_full <- as.data.frame(fread(file.path(input_dir, "genome_wide_marker_map.tsv")))
GM <- GM_full[, c("SNP", "Chromosome", "Position")]

## ---- verbatim from 01_fit_m4_exact_partitions.R ----
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
## ---- end verbatim ----

y <- Y[[trait]]
mask <- fold_mask(fold); y_train <- y; y_train[c(mask$test,mask$excluded)] <- NA
train_rows <- which(!is.na(y_train))
train_df <- Y[train_rows,c("strain","env",trait)]
envs <- sort(unique(train_df$env))
env_results <- lapply(envs, function(e) run_env_gwas(train_df,e,fold)); names(env_results) <- envs
ranking <- as.data.table(meta_rank(env_results)$ranking)
fwrite(ranking, file.path(out_dir, "meta_analysis_ranking.csv.gz"))

## ---- plot, as in 03_gwas_plots_and_tables.R ----
d <- merge(ranking, as.data.table(GM_full), by="SNP", all.x=TRUE, sort=FALSE)
setorder(d,P_meta,-N_environments,min_P,SNP)
d[,FDR_meta:=p.adjust(P_meta,method="BH")]
d[,Rank:=seq_len(.N)]
draw <- function() {
  chr_order <- paste0("Ca",1:8)
  pdat <- d[Chromosome %in% chr_order & is.finite(Position) & is.finite(P_meta)]
  pdat[,Chromosome:=factor(Chromosome,levels=chr_order)]
  chr_len <- pdat[,.(len=max(Position,na.rm=TRUE)),keyby=Chromosome]  # keyby: lay out Ca1..Ca8 in order
  chr_len[,offset:=shift(cumsum(len),fill=0)]
  pdat <- merge(pdat,chr_len[,.(Chromosome,offset)],by="Chromosome",sort=FALSE)
  pdat[,x:=Position+offset]
  cols <- rep(c("#2C7FB8","#7FCDBB"),4)
  plot(pdat$x,-log10(pmax(pdat$P_meta,.Machine$double.xmin)),pch=20,cex=.28,
       col=cols[as.integer(pdat$Chromosome)],xaxt="n",xlab="Chromosome",ylab=expression(-log[10](italic(P))),
       main=sprintf("%s %s partition %d fold %d: meta-GWAS",trait,cv_name,partition,fold))
  axis(1,at=pdat[,mean(range(x)),keyby=Chromosome]$V1,labels=chr_order)
  bonf_p <- 0.05 / sum(is.finite(d$P_meta))
  abline(h=-log10(bonf_p),col="#D7301F",lty=2,lwd=1.5)
  bh_pass <- d[is.finite(P_meta) & FDR_meta <= 0.05, P_meta]
  if (length(bh_pass)) abline(h=-log10(max(bh_pass)),col="#762A83",lty=3,lwd=1.5)
  legend("topright",
         legend=c("Bonferroni 0.05",if(length(bh_pass)) "BH-FDR 5%"),
         col=c("#D7301F",if(length(bh_pass)) "#762A83"),
         lty=c(2,if(length(bh_pass)) 3),lwd=1.5,bty="n",cex=.8)
  obs <- sort(pmax(pmin(d$P_meta,1),.Machine$double.xmin))
  exp <- ppoints(length(obs))
  plot(-log10(exp),-log10(obs),pch=20,cex=.35,xlab=expression(Expected~~-log[10](italic(P))),
       ylab=expression(Observed~~-log[10](italic(P))),
       main=sprintf("%s %s partition %d fold %d: QQ",trait,cv_name,partition,fold))
  abline(0,1,col="red",lwd=1.5)
}
stem <- file.path(out_dir, sprintf("%s_%s_partition%02d_fold%d", trait, cv_name, partition, fold))
pdf(paste0(stem, ".pdf"), width=12, height=6, onefile=TRUE); draw(); dev.off()
png(paste0(stem, "_%d.png"), width=12, height=6, units="in", res=150); draw(); dev.off()

## ---- reproducibility check against the deposited top 500 ----
ref <- fread(file.path(project_dir, "reference_results", "gwas", "selected_markers_all_folds.tsv"))
ref <- ref[Partition == partition & Trait == trait & CV == cv_name & Fold == fold]
new <- d[Rank <= 500]
chk <- merge(ref[, .(SNP, Rank_ref = Rank, P_ref = P_meta, Z_ref = Z_meta)],
             new[, .(SNP, Rank_new = Rank, P_new = P_meta, Z_new = Z_meta)], by = "SNP", all = TRUE)
res <- data.table(partition, trait, CV = cv_name, fold, n_env = length(envs),
                  overlap = chk[!is.na(Rank_ref) & !is.na(Rank_new), .N],
                  same_order = identical(ref[order(Rank), SNP], new[order(Rank), SNP]),
                  max_rel_diff_Z = chk[!is.na(Z_ref) & !is.na(Z_new), max(abs(Z_new - Z_ref) / abs(Z_ref))])
fwrite(res, file.path(out_dir, "reproducibility_check.csv"))
print(res)
file.create(done)
