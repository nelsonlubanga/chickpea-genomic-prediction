#!/usr/bin/env Rscript

## Redraw the Manhattan/Q-Q plots of one fold from its saved ranking
## (meta_analysis_ranking.csv.gz in the fold folder),
## without rerunning the GWAS. Overwrites the fold's .pdf and _1/_2.png.
## Usage: Rscript 11_redraw_fold_gwas_plots.R <fold_dir> [<marker_map.tsv>]

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
fold_dir <- normalizePath(args[1], mustWork = TRUE)
map_file <- if (length(args) >= 2L) args[2] else file.path(Sys.getenv("PROJECT_DIR", unset = "."), "inputs", "genome_wide_marker_map.tsv")

parts <- rev(strsplit(fold_dir, "/", fixed = TRUE)[[1]])[1:4]
fold <- as.integer(sub("fold_", "", parts[1])); cv_name <- parts[2]; trait <- parts[3]
partition <- as.integer(sub("partition_", "", parts[4]))

ranking <- fread(file.path(fold_dir, "meta_analysis_ranking.csv.gz"))
GM_full <- fread(map_file)

## ---- identical to the plot block of 10_regenerate_fold_gwas_plots.R ----
d <- merge(ranking, GM_full, by="SNP", all.x=TRUE, sort=FALSE)
setorder(d,P_meta,-N_environments,min_P,SNP)
d[,FDR_meta:=p.adjust(P_meta,method="BH")]
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
stem <- file.path(fold_dir, sprintf("%s_%s_partition%02d_fold%d", trait, cv_name, partition, fold))
## REDRAW_PDF=0 writes the PNGs only (each vector PDF is ~20 MB).
if (Sys.getenv("REDRAW_PDF", unset = "1") != "0") {
  pdf(paste0(stem, ".pdf"), width=12, height=6, onefile=TRUE); draw(); dev.off()
}
png(paste0(stem, "_%d.png"), width=12, height=6, units="in", res=150); draw(); dev.off()
