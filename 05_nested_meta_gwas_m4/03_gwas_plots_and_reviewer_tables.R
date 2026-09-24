#!/usr/bin/env Rscript

## Post-process completed Alemu-style folds into reviewer-facing GWAS evidence.
## Run only after the Slurm prediction array has completed.

suppressPackageStartupMessages(library(data.table))
project <- normalizePath(Sys.getenv("PROJECT_DIR", unset="."), mustWork=TRUE)
results <- file.path(project,"results")
out <- file.path(project,"reviewer_outputs")
dir.create(file.path(out,"gwas_plots"),recursive=TRUE,showWarnings=FALSE)
map <- fread(file.path(project,"inputs","genome_wide_marker_map.tsv"))

ranking_files <- list.files(results,pattern="meta_analysis_ranking\\.csv$",recursive=TRUE,full.names=TRUE)
if (!length(ranking_files)) stop("No completed meta-analysis rankings found")

parse_key <- function(f) {
  x <- strsplit(sub(paste0("^",results,"/?"),"",f),"/",fixed=TRUE)[[1]]
  data.table(Partition=as.integer(sub("partition_","",x[1])),Trait=x[2],CV=x[3],Fold=as.integer(sub("fold_","",x[4])))
}

selected_all <- vector("list",length(ranking_files))
for (i in seq_along(ranking_files)) {
  f <- ranking_files[i]; key <- parse_key(f); d <- fread(f)
  d <- merge(d,map,by="SNP",all.x=TRUE,sort=FALSE)
  setorder(d,P_meta,-N_environments,min_P,SNP)
  d[,FDR_meta:=p.adjust(P_meta,method="BH")]
  d[,Rank:=seq_len(.N)]
  selected_all[[i]] <- cbind(key[rep(1,min(500,.N))],d[seq_len(min(500,.N))])

  pdf_name <- sprintf("%s_%s_partition%02d_fold%d.pdf",key$Trait,key$CV,key$Partition,key$Fold)
  pdf(file.path(out,"gwas_plots",pdf_name),width=12,height=6,onefile=TRUE)
  chr_order <- paste0("Ca",1:8)
  pdat <- d[Chromosome %in% chr_order & is.finite(Position) & is.finite(P_meta)]
  pdat[,Chromosome:=factor(Chromosome,levels=chr_order)]
  chr_len <- pdat[,.(len=max(Position,na.rm=TRUE)),by=Chromosome]
  chr_len[,offset:=shift(cumsum(len),fill=0)]
  pdat <- merge(pdat,chr_len[,.(Chromosome,offset)],by="Chromosome",sort=FALSE)
  pdat[,x:=Position+offset]
  cols <- rep(c("#2C7FB8","#7FCDBB"),4)
  plot(pdat$x,-log10(pmax(pdat$P_meta,.Machine$double.xmin)),pch=20,cex=.28,
       col=cols[as.integer(pdat$Chromosome)],xaxt="n",xlab="Chromosome",ylab=expression(-log[10](italic(P))),
       main=sprintf("%s %s partition %d fold %d: meta-GWAS",key$Trait,key$CV,key$Partition,key$Fold))
  axis(1,at=pdat[,mean(range(x)),by=Chromosome]$V1,labels=chr_order)
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
       main=sprintf("%s %s partition %d fold %d: QQ",key$Trait,key$CV,key$Partition,key$Fold))
  abline(0,1,col="red",lwd=1.5)
  dev.off()
}

selected <- rbindlist(selected_all,fill=TRUE)
fwrite(selected,file.path(out,"selected_markers_all_folds.tsv"),sep="\t")
freq <- selected[,.(Selections=.N,Partitions=uniqueN(Partition),Folds=uniqueN(paste(Partition,Fold))),
                 by=.(Trait,CV,SNP,Chromosome,Position,A1,A2)]
freq[,Selection_frequency:=Selections/50]
setorder(freq,Trait,CV,-Selections,SNP)
fwrite(freq,file.path(out,"marker_selection_frequency_across_partitions.tsv"),sep="\t")

env_files <- list.files(results,pattern="selected_marker_environment_results\\.csv$",recursive=TRUE,full.names=TRUE)
env_all <- rbindlist(lapply(env_files,function(f)cbind(parse_key(f)[rep(1,nrow(fread(f)))],fread(f))),fill=TRUE)
env_all <- merge(env_all,map,by="SNP",all.x=TRUE,sort=FALSE)
fwrite(env_all,file.path(out,"selected_marker_environment_evidence.tsv"),sep="\t")

audit <- data.table(
  Item=c("Marker discovery","Population structure","Leakage prevention","Across-environment evidence",
         "Marker selection","Selected-marker model","Proxy rule","GWAS plots","Selection frequency"),
  Specification=c("BLINK in GAPIT 4.1.0","First three marker principal components",
    "GWAS repeated within every training fold; validation phenotypes excluded",
    "Signed sample-size-weighted Stouffer meta-analysis","Top 500 ranked meta-analysis SNPs",
    "BRR-shrunk random effects in E+L+G+GxE+S","Best tag within 500 kb only if r2 >= 0.5; duplicates collapsed",
    "Manhattan and QQ PDF for every trait x CV x partition x fold","Reported across 10 partitions and 5 folds"))
fwrite(audit,file.path(out,"GWAS_reproducibility_audit.tsv"),sep="\t")

pred_files <- list.files(results,pattern="^M4_(chr|tags)_top500\\.csv$",recursive=TRUE,full.names=TRUE)
completeness <- rbindlist(lapply(pred_files,function(f) {
  key <- parse_key(f); z <- fread(f)
  panel <- sub("^M4_(chr|tags)_top500\\.csv$","\\1",basename(f))
  cbind(key,Panel=panel,N_predictions=nrow(z),N_environments=uniqueN(z$environment),
        Missing_y=sum(is.na(z$y)),Missing_yHat=sum(is.na(z$yHat)))
}),fill=TRUE)
fwrite(completeness,file.path(out,"M4_cross_validation_completeness.tsv"),sep="\t")
expected <- CJ(Partition=1:10,Trait=c("DTF","DTM","HSW","PH","PPP","YPPlnt"),
               CV=c("CV0","CV1","CV2","CV00"),Fold=1:5,Panel=c("chr","tags"))
missing <- fsetdiff(expected,unique(completeness[,.(Partition,Trait,CV,Fold,Panel)]))
fwrite(missing,file.path(out,"M4_missing_expected_results.tsv"),sep="\t")
cat("Created reviewer outputs in",out,"\n")
