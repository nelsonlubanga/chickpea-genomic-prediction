#!/usr/bin/env Rscript

## Sensitivity analysis reusing completed training-only meta-GWAS rankings.
## Usage: Rscript 08_fit_top100_top1000.R <partition> <trait> <CV>

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: Rscript 08_fit_top100_top1000.R <partition> <trait> <CV>")
partition <- as.integer(args[1]); trait <- args[2]; cv_name <- args[3]
traits <- c("DTF","DTM","HSW","PH","PPP","YPPlnt")
cvs <- c("CV0","CV1","CV2","CV00")
stopifnot(partition%in%1:10,trait%in%traits,cv_name%in%cvs)
suppressPackageStartupMessages({library(BGLR);library(data.table)})

project <- normalizePath(Sys.getenv("PROJECT_DIR",unset="."),mustWork=TRUE)
input <- file.path(project,"inputs"); results <- file.path(project,"results")
scratch <- Sys.getenv("SLURM_TMPDIR",unset="")
if (!nzchar(scratch)) scratch <- file.path(tempdir(),"m4_sensitivity_tmp")
dir.create(scratch,recursive=TRUE,showWarnings=FALSE)
n_iter <- 5000L; burn_in <- 2000L; thin <- 5L
ks <- c(100L,1000L); ld_window <- 500000; min_r2 <- 0.5
set.seed(22345+partition)

obj <- readRDS(file.path(input,"step5_panels.rds"))
Y <- obj$Y; E <- obj$E; L <- obj$L; panels <- obj$panels
accessions <- readRDS(file.path(input,"step5_base.rds"))$accessions
p12 <- readRDS(file.path(input,"step5_partitions.rds"))[[partition]]
p00 <- readRDS(file.path(input,"step5_partitions_cv0_cv00.rds"))[[partition]]
GD <- readRDS(file.path(input,"genome_wide_genotypes_190x409128.rds"))
GD <- GD[match(accessions,GD$Taxa),,drop=FALSE]
X <- as.matrix(GD[-1]);rownames(X)<-accessions;storage.mode(X)<-"double"
for(j in which(colSums(is.na(X))>0L)) X[is.na(X[,j]),j] <- mean(X[,j],na.rm=TRUE)
GM <- as.data.frame(fread(file.path(input,"genome_wide_marker_map.tsv")));rownames(GM)<-GM$SNP
tags <- scan(file.path(input,"haplotype_tagged_marker_list.txt"),what="",quiet=TRUE)

fold_mask <- function(fold) {
  if(cv_name=="CV1") return(list(test=which(p12$CV1==fold),excluded=integer()))
  if(cv_name=="CV2") return(list(test=which(p12$CV2==fold),excluded=integer()))
  if(cv_name=="CV0") return(list(test=which(p00$CV0==fold),excluded=integer()))
  list(test=which(p00$geno_fold==fold&p00$env_fold==fold),
       excluded=which(xor(p00$geno_fold==fold,p00$env_fold==fold)))
}

map_ranked_to_tags <- function(snps) {
  out <- vector("list",length(snps))
  for(i in seq_along(snps)) {
    s <- snps[i]
    if(s%in%tags) {out[[i]]<-data.frame(rank=i,discovery=s,proxy=s,r2=1);next}
    if(!s%in%rownames(GM)) next
    chr<-as.character(GM[s,"Chromosome"]);pos<-as.numeric(GM[s,"Position"])
    cand<-tags[startsWith(tags,paste0(chr,"_"))]
    if(!length(cand)) next
    cp<-suppressWarnings(as.numeric(sub("^.*_","",cand)))
    near<-cand[is.finite(cp)&abs(cp-pos)<=ld_window]
    if(!length(near)) next
    r2<-vapply(near,function(tg)suppressWarnings(cor(X[,s],X[,tg])^2),numeric(1))
    r2[!is.finite(r2)]<-0;b<-which.max(r2)
    if(r2[b]>=min_r2) out[[i]]<-data.frame(rank=i,discovery=s,proxy=near[b],r2=r2[b])
  }
  ans<-rbindlist(out,fill=TRUE)
  if(!nrow(ans)) ans<-data.table(rank=integer(),discovery=character(),proxy=character(),r2=numeric())
  ans
}

fit_model <- function(markers,panel_name,k,y_train,y,test,out_dir) {
  p<-panels[[panel_name]]
  ETA<-list(list(X=E,model="BRR"),list(X=L,model="BRR"),
            list(V=p$EVD.G$vectors,d=p$EVD.G$values,model="RKHS"),
            list(V=p$EVD.GxE$vectors,d=p$EVD.GxE$values,model="RKHS"))
  if(length(markers)) {
    xs<-X[,markers,drop=FALSE];xs<-xs[,apply(xs,2,sd)>0,drop=FALSE]
    if(ncol(xs)) ETA[[5]]<-list(X=scale(L%*%xs),model="BRR")
  }
  prefix<-file.path(scratch,sprintf("sens_p%d_%s_%s_%s_k%d_",partition,trait,cv_name,panel_name,k))
  fm<-BGLR(y=y_train,ETA=ETA,nIter=n_iter,burnIn=burn_in,thin=thin,verbose=FALSE,saveAt=prefix)
  unlink(Sys.glob(paste0(prefix,"*.dat")))
  test<-test[!is.na(y[test])]
  fwrite(data.table(testing=test,Individual=Y$strain[test],environment=Y$env[test],
                    y=y[test],yHat=fm$yHat[test]),
         file.path(out_dir,sprintf("M4_%s_top%d.csv",panel_name,k)))
}

y<-Y[[trait]]
for(fold in 1:5) {
  out_dir<-file.path(results,sprintf("partition_%02d",partition),trait,cv_name,paste0("fold_",fold))
  done<-file.path(out_dir,".sensitivity_top100_top1000.complete")
  if(file.exists(done)){message("Skipping completed ",out_dir);next}
  ranking_file<-file.path(out_dir,"meta_analysis_ranking.csv")
  if(!file.exists(ranking_file))stop("Missing ranking: ",ranking_file)
  ranking<-fread(ranking_file);setorder(ranking,P_meta,-N_environments,min_P,SNP)
  top1000<-head(ranking$SNP,1000L)
  tag_map<-map_ranked_to_tags(top1000)
  fwrite(tag_map,file.path(out_dir,"tag_proxy_mapping_top1000.csv"))
  mask<-fold_mask(fold);y_train<-y;y_train[c(mask$test,mask$excluded)]<-NA
  for(k in ks) {
    selected<-head(top1000,k);writeLines(selected,file.path(out_dir,sprintf("selected_top%d.txt",k)))
    fit_model(intersect(selected,colnames(X)),"chr",k,y_train,y,mask$test,out_dir)
    proxies<-if(nrow(tag_map)) unique(tag_map[rank<=k]$proxy) else character()
    fit_model(proxies,"tags",k,y_train,y,mask$test,out_dir)
    if(!length(proxies)) writeLines(sprintf("No qualifying tag proxy among top %d; M4 reduced to M3.",k),
                                    file.path(out_dir,sprintf("M4_tags_top%d_M3_fallback.txt",k)))
  }
  file.create(done)
  message(sprintf("p%d %s %s fold%d sensitivity complete",partition,trait,cv_name,fold))
}
