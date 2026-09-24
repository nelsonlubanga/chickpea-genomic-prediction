project_dir <- path.expand(Sys.getenv("PROJECT_DIR", unset = "."))
partition_dir <- path.expand(Sys.getenv("PARTITION_DIR", unset = "exact_partitions"))
setwd(project_dir)
traits  <- c("DTF", "DTM", "HSW", "PH", "PPP", "YPPlnt")
n_folds <- 5

p_cv12   <- readRDS(file.path(partition_dir, "step5_partitions.rds"))
p_cv0000 <- readRDS(file.path(partition_dir, "step5_partitions_cv0_cv00.rds"))

n_checked <- 0; n_fail <- 0; total_extra <- 0

check_subset <- function(p, cv_name, fold, trait, recon_testing) {
  f <- file.path("output_mechanism", paste0("p", p), trait, cv_name, paste0("fold_", fold), "M1.csv")
  if (!file.exists(f)) return(invisible(NULL))
  d <- read.csv(f)
  ok <- all(d$testing %in% recon_testing)
  n_checked <<- n_checked + 1
  if (!ok) {
    n_fail <<- n_fail + 1
    cat(sprintf("  FAIL p%d %s %s fold%d: %d/%d ground-truth testing rows NOT in reconstruction\n",
                p, trait, cv_name, fold, sum(!(d$testing %in% recon_testing)), length(d$testing)))
  }
}

for (p in 1:10) {
  cv12   <- p_cv12[[p]]
  cv0000 <- p_cv0000[[p]]
  for (trait in traits) {
    for (fold in 1:n_folds) {
      check_subset(p, "CV1", fold, trait, which(cv12$CV1 == fold))
      check_subset(p, "CV2", fold, trait, which(cv12$CV2 == fold))
      check_subset(p, "CV0", fold, trait, which(cv0000$CV0 == fold))
      ## raw CSV `testing` column for CV00 = true_test UNION excluded =
      ## {geno_fold==fold} UNION {env_fold==fold} -- see reconstruct_step5.R
      check_subset(p, "CV00", fold, trait, which(cv0000$geno_fold == fold | cv0000$env_fold == fold))
    }
  }
}

cat(sprintf("\nChecked %d (partition x trait x cv x fold) ground-truth files. Failures: %d\n", n_checked, n_fail))
if (n_checked == 0) {
  stop("No ground-truth M1.csv files were found under '", partition_dir,
       "/output_mechanism' -- this check validates NOTHING until M1-M3 has ",
       "actually been fit and its output placed there. A silent 0-checked, ",
       "0-failed result must not be reported as a pass.")
}
if (n_fail == 0) cat("ALL PASS: every already-completed M1.csv testing set is a subset of the reconstructed fold. Reconstruction is consistent with all 10 partitions' existing M1-M3 output.\n")
