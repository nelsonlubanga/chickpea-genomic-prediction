#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
project <- normalizePath(Sys.getenv("PROJECT_DIR",unset="."),mustWork=TRUE)
results <- file.path(project,"results")
out <- file.path(project,"reviewer_outputs")
parse_key <- function(f) {
  x <- strsplit(sub(paste0("^",results,"/?"),"",f),"/",fixed=TRUE)[[1]]
  data.table(Partition=as.integer(sub("partition_","",x[1])),Trait=x[2],CV=x[3],Fold=as.integer(sub("fold_","",x[4])))
}
fs <- list.files(results,pattern="^M4_(chr|tags)_top500\\.csv$",recursive=TRUE,full.names=TRUE)
z <- rbindlist(lapply(fs,function(f) {
  d <- fread(f)
  panel <- sub("^M4_(chr|tags)_top500\\.csv$","\\1",basename(f))
  cbind(parse_key(f),Panel=panel,N_predictions=nrow(d),N_environments=uniqueN(d$environment),
        Missing_y=sum(is.na(d$y)),Missing_yHat=sum(is.na(d$yHat)))
}))
fwrite(z,file.path(out,"M4_cross_validation_completeness.tsv"),sep="\t")
expected <- CJ(Partition=1:10,Trait=c("DTF","DTM","HSW","PH","PPP","YPPlnt"),
               CV=c("CV0","CV1","CV2","CV00"),Fold=1:5,Panel=c("chr","tags"))
missing <- fsetdiff(expected,unique(z[,.(Partition,Trait,CV,Fold,Panel)]))
fwrite(missing,file.path(out,"M4_missing_expected_results.tsv"),sep="\t")
cat("Prediction files:",nrow(z),"; missing expected files:",nrow(missing),"\n")
