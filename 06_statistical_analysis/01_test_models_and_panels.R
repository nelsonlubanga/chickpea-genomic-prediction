suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(emmeans); library(writexl)
})

args <- commandArgs(trailingOnly=TRUE)
top_n <- if (length(args) >= 1L) as.integer(args[[1]]) else 500L
infile <- if (length(args) >= 2L) args[[2L]] else
  sprintf("reference_results/data/PA_fold_level_M1_M4_top%d.csv", top_n)
outfile <- if (length(args) >= 3L) args[[3L]] else
  sprintf("outputs/Table_M4_top%d_model_panel_comparison.xlsx", top_n)
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
if (!file.exists(infile)) stop("Input file not found: ", infile)
d <- read.csv(infile)
d$Model <- factor(d$Model, levels=c("M1","M2","M3","M4"))
d$Panel <- factor(d$Panel, levels=c("genome-wide","haplotype-tagged"))
d$Repetition <- factor(d$Repetition); d$Fold <- factor(d$Fold)

traits <- c("DTF","DTM","HSW","PH","PPP","YPPlnt")
cvs <- c("CV1","CV2","CV0","CV00")
panels <- levels(d$Panel)

mc_o <- list(); mc_p <- list(); pc_o <- list(); pc_p <- list()
for (tr in traits) for (cv in cvs) for (pn in panels) {
  z <- droplevels(subset(d, Trait==tr & CV==cv & Panel==pn & is.finite(PA)))
  fit <- suppressWarnings(lmer(PA~Model+(1|Repetition/Fold),z,REML=TRUE,
                               control=lmerControl(optimizer="bobyqa")))
  a <- anova(fit)
  mc_o[[length(mc_o)+1]] <- data.frame(Trait=tr,CV=cv,Panel=pn,F=a$`F value`[1],
    NumDF=a$NumDF[1],DenDF=a$DenDF[1],P=a$`Pr(>F)`[1])
  if (a$`Pr(>F)`[1] < .05) {
    x <- as.data.frame(pairs(emmeans(fit,~Model),adjust="tukey"))
    x$Trait<-tr;x$CV<-cv;x$Panel<-pn;mc_p[[length(mc_p)+1]]<-x
  }
}
for (tr in traits) for (cv in cvs) for (mo in c("M2","M3","M4")) {
  z <- droplevels(subset(d,Trait==tr & CV==cv & Model==mo & is.finite(PA)))
  fit <- suppressWarnings(lmer(PA~Panel+(1|Repetition/Fold),z,REML=TRUE,
                               control=lmerControl(optimizer="bobyqa")))
  a <- anova(fit); x <- as.data.frame(pairs(emmeans(fit,~Panel),adjust="none"))
  pc_o[[length(pc_o)+1]] <- data.frame(Trait=tr,CV=cv,Model=mo,F=a$`F value`[1],
    NumDF=a$NumDF[1],DenDF=a$DenDF[1],P=a$`Pr(>F)`[1])
  x$Trait<-tr;x$CV<-cv;x$Model<-mo;pc_p[[length(pc_p)+1]]<-x
}
mc_o<-do.call(rbind,mc_o);mc_p<-do.call(rbind,mc_p)
pc_o<-do.call(rbind,pc_o);pc_p<-do.call(rbind,pc_p)

fix_trait <- function(x) {x[x=="YPPlnt"]<-"YPP";x}
fix_panel <- function(x) ifelse(x=="genome-wide","Genome-wide","Haplotype-tagged")
fix_model <- function(x) {x<-as.character(x);x[x=="M4"]<-sprintf("M4-%d",top_n);x}
ord <- function(x,extra) {
  x$Trait<-factor(x$Trait,levels=c("DTF","DTM","HSW","PH","PPP","YPP"))
  x$CV<-factor(x$CV,levels=c("CV1","CV2","CV0","CV00"))
  x[order(x$Trait,x$CV,x[[extra]]),]
}

mc_om <- data.frame(Trait=fix_trait(mc_o$Trait),CV=mc_o$CV,Panel=fix_panel(mc_o$Panel),
  Comparison=sprintf("M1 vs M2 vs M3 vs M4-%d",top_n),F=round(mc_o$F,2),`Num DF`=mc_o$NumDF,
  `Den DF`=round(mc_o$DenDF,1),`p-value`=mc_o$P,check.names=FALSE)
mc_om<-ord(mc_om,"Panel")
mc_pw <- data.frame(Trait=fix_trait(mc_p$Trait),CV=mc_p$CV,Panel=fix_panel(mc_p$Panel),
  Contrast=gsub("M4",sprintf("M4-%d",top_n),mc_p$contrast,fixed=TRUE),
  `Estimate (PA difference)`=round(mc_p$estimate,4),SE=round(mc_p$SE,4),
  DF=mc_p$df,`t-ratio`=round(mc_p$t.ratio,4),
  `p-value (Tukey-adjusted)`=mc_p$p.value,check.names=FALSE)
mc_pw<-ord(mc_pw,"Panel")
pc_om <- data.frame(Trait=fix_trait(pc_o$Trait),CV=pc_o$CV,Model=fix_model(pc_o$Model),
  Comparison="Genome-wide vs haplotype-tagged",F=round(pc_o$F,2),
  `Num DF`=pc_o$NumDF,`Den DF`=round(pc_o$DenDF,1),`p-value`=pc_o$P,check.names=FALSE)
pc_om<-ord(pc_om,"Model")
pc_pw <- data.frame(Trait=fix_trait(pc_p$Trait),CV=pc_p$CV,Model=fix_model(pc_p$Model),
  Contrast="Genome-wide - Haplotype-tagged",
  `Estimate (PA difference)`=round(pc_p$estimate,4),SE=round(pc_p$SE,4),DF=pc_p$df,
  `t-ratio`=round(pc_p$t.ratio,4),`p-value`=pc_p$p.value,check.names=FALSE)
pc_pw<-ord(pc_pw,"Model")

readme <- data.frame(
  col1=c(sprintf("Top-%d sensitivity analysis: statistical comparison of predictive ability among M1, M2, M3 and M4-%d and between marker panels, by trait and cross-validation scheme.",top_n,top_n),NA,
         "Sheet","Model_comparison_omnibus","Model_comparison_pairwise","Panel_comparison_omnibus","Panel_comparison_pairwise",NA,
         "Statistical model","M4 specification","Predictive ability","Pairwise contrasts"),
  col2=c(NA,NA,"Contents",sprintf("Omnibus F-test among M1, M2, M3 and M4-%d per trait x CV x panel.",top_n),
         "Tukey-adjusted pairwise contrasts, reported where the omnibus model test was significant (p < 0.05).",
         sprintf("Omnibus genome-wide versus haplotype-tagged test within M2, M3 and M4-%d.",top_n),
         "Genome-wide minus haplotype-tagged contrast for every trait x CV x model combination.",NA,
         "Fold-level PA ~ Model or Panel + (1 | Repetition/Fold), fitted separately by trait and CV scheme.",
         sprintf("M4-%d contains the %d highest-ranked fold-specific GWAS markers; M1-M3 are unchanged from the primary analysis.",top_n,top_n),
         "Within-environment Pearson correlations for environments with at least five validation observations, averaged equally within fold.",
         "emmeans; Tukey adjustment for six model contrasts and no additional adjustment for the single panel contrast."),
  check.names=FALSE)
names(readme)<-c(sprintf("Top-%d marker-set sensitivity analysis",top_n),"")

write_xlsx(list(README=readme,Model_comparison_omnibus=mc_om,
  Model_comparison_pairwise=mc_pw,Panel_comparison_omnibus=pc_om,
  Panel_comparison_pairwise=pc_pw),outfile)
cat("Wrote",outfile,"\n",nrow(mc_om),nrow(mc_pw),nrow(pc_om),nrow(pc_pw),"rows\n")
