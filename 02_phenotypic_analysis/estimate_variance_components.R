library(asreml)

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else "inputs/second_stage_BLUEs_Y.csv"
output_file <- if (length(args) >= 2L) args[[2L]] else
  "outputs/variance_components_heritability.csv"
## Optional third argument: file for the full REML variance-component table.
varcomp_file <- if (length(args) >= 3L) args[[3L]] else NA_character_
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dat <- read.csv(input_file, stringsAsFactors = FALSE)
dat <- transform(dat, Genotype = factor(strain), year = factor(pwdyear),
                 location = factor(location), env = factor(env))
traits <- c("DTF", "PH", "DTM", "PPP", "HSW", "YPPlnt")
out <- list()
full <- list()
for (tr in traits) {
  d <- dat[, c("Genotype", "env", "year", "location", tr)]
  names(d)[5] <- "y"
  d <- d[is.finite(d$y), ]
  fit <- asreml(fixed = y ~ 1 + env,
               random = ~ Genotype + Genotype:location + Genotype:year,
               residual = ~ id(units), data = d, trace = FALSE)
  for (i in seq_len(10)) {
    if (fit$converge) break
    fit <- update(fit)
  }
  vc <- summary(fit)$varcomp
  getvc <- function(pattern) {
    z <- vc[grepl(pattern, rownames(vc), fixed = FALSE), "component"]
    if (length(z)) as.numeric(z[1]) else NA_real_
  }
  vg <- getvc("^Genotype$")
  vgs <- getvc("Genotype:location")
  vgy <- getvc("Genotype:year")
  ve <- getvc("units!")
  nobs <- mean(table(d$Genotype)); s <- nlevels(d$location); y <- nlevels(d$year)
  h2 <- vg / (vg + vgs / s + vgy / y + ve / nobs)
  out[[tr]] <- data.frame(trait = tr, converged = fit$converge,
    genetic_variance = vg, Gxsite_variance = vgs, Gxyear_variance = vgy,
    residual_variance = ve, mean_observations_per_genotype = nobs,
    sites = s, years = y, H2_entry_mean = h2)
  full[[tr]] <- data.frame(trait = tr, term = rownames(vc),
    vc[, c("component", "std.error", "z.ratio", "bound")],
    pct_of_total = 100 * vc$component / sum(vc$component),
    n_obs = nrow(d), n_lines = nlevels(droplevels(d$Genotype)),
    logLik = fit$loglik, converged = fit$converge, row.names = NULL)
}
res <- do.call(rbind, out)
write.csv(res, output_file, row.names = FALSE)
if (!is.na(varcomp_file)) {
  dir.create(dirname(varcomp_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(do.call(rbind, full), varcomp_file, row.names = FALSE)
}
print(res)
