suppressPackageStartupMessages(library(BGLR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop(paste(
  "Usage: Rscript 01_fit_m1_m3_exact_partitions.R",
  "<kernel_objects.rds> <cv12.rds> <cv0_cv00.rds> <partition 1..10> <output_dir>"))
obj <- readRDS(args[1]); cv12_all <- readRDS(args[2]); cv00_all <- readRDS(args[3])
partition <- as.integer(args[4]); out_root <- args[5]
stopifnot(partition %in% 1:10)
Y <- obj$Y; E <- obj$E; L <- obj$L; panels <- obj$panels
cv12 <- cv12_all[[partition]]; cv00 <- cv00_all[[partition]]
available_traits <- c("DTF", "DTM", "HSW", "PH", "PPP", "YPPlnt")
traits_arg <- Sys.getenv("TRAITS", unset = paste(available_traits, collapse = ","))
traits <- trimws(strsplit(traits_arg, ",", fixed = TRUE)[[1]])
if (!length(traits) || any(!traits %in% available_traits)) {
  stop("TRAITS must be a comma-separated subset of: ",
       paste(available_traits, collapse = ", "))
}
available_cv <- c("CV1", "CV2", "CV0", "CV00")
cv_arg <- Sys.getenv("CV_SCHEMES", unset = paste(available_cv, collapse = ","))
cv_schemes <- trimws(strsplit(cv_arg, ",", fixed = TRUE)[[1]])
if (!length(cv_schemes) || any(!cv_schemes %in% available_cv)) {
  stop("CV_SCHEMES must be a comma-separated subset of: ",
       paste(available_cv, collapse = ", "))
}
set.seed(12345)

mask <- function(cv, fold) {
  if (cv == "CV1") return(list(test = which(cv12$CV1 == fold), excluded = integer()))
  if (cv == "CV2") return(list(test = which(cv12$CV2 == fold), excluded = integer()))
  if (cv == "CV0") return(list(test = which(cv00$CV0 == fold), excluded = integer()))
  test <- which(cv00$geno_fold == fold & cv00$env_fold == fold)
  excluded <- which(xor(cv00$geno_fold == fold, cv00$env_fold == fold))
  list(test = test, excluded = excluded)
}

bglr_scratch <- file.path(tempdir(), "m1_m3_bglr")
dir.create(bglr_scratch, recursive = TRUE, showWarnings = FALSE)

fit_one <- function(y, m, panel = NULL, tag = "run") {
  eta <- list(list(X = E, model = "BRR"), list(X = L, model = "BRR"))
  if (m %in% c("M2", "M3")) {
    p <- panels[[panel]]
    eta <- c(eta, list(list(V = p$EVD.G$vectors, d = p$EVD.G$values, model = "RKHS")))
    if (m == "M3") eta <- c(eta, list(list(V = p$EVD.GxE$vectors, d = p$EVD.GxE$values, model = "RKHS")))
  }
  ## Every call needs its own saveAt prefix -- without one, BGLR writes its
  ## trace files (mu.dat, varE.dat, ETA_*.dat) to the current working
  ## directory using the same fixed names every time, so repeated calls
  ## (every trait x CV x fold x model x panel combination) silently
  ## overwrite each other and litter wherever the script happens to be run
  ## from.
  prefix <- file.path(bglr_scratch, paste0(tag, "_"))
  fm <- BGLR(y = y, ETA = eta, nIter = 5000, burnIn = 2000, thin = 5,
             verbose = FALSE, saveAt = prefix)
  unlink(Sys.glob(paste0(prefix, "*.dat")))
  fm
}

for (trait in traits) for (cv in cv_schemes) for (fold in 1:5) {
  y <- Y[[trait]]; z <- mask(cv, fold); y_train <- y
  y_train[c(z$test, z$excluded)] <- NA
  reported_test <- z$test[!is.na(y[z$test])]
  out <- file.path(out_root, paste0("p", partition), trait, cv, paste0("fold_", fold))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  run_tag <- paste("p", partition, trait, cv, fold, sep = "_")
  fm <- fit_one(y_train, "M1", tag = paste0(run_tag, "_M1"))
  write.csv(data.frame(testing = reported_test, Individual = Y$strain[reported_test],
                       environment = Y$env[reported_test],
                       y = y[reported_test], yHat = fm$yHat[reported_test]),
            file.path(out, "M1.csv"), row.names = FALSE)
  for (model in c("M2", "M3")) for (panel in c("chr", "tags")) {
    fm <- fit_one(y_train, model, panel, tag = paste0(run_tag, "_", model, "_", panel))
    write.csv(data.frame(testing = reported_test, Individual = Y$strain[reported_test],
                         environment = Y$env[reported_test],
                         y = y[reported_test], yHat = fm$yHat[reported_test]),
              file.path(out, paste0(model, "_", panel, ".csv")), row.names = FALSE)
  }
}
