## First-stage analysis described in the manuscript. This script documents
## the model ICRISAT used to obtain within-environment BLUEs from their own
## raw plot-level field records (which include block assignments this team
## has never had access to) -- it is NOT executed by this pipeline. The
## resulting BLUEs were shared with this study directly as
## inputs/second_stage_BLUEs_Y.csv; the raw plot_level_phenotypes.csv this
## script expects does not exist in this team's data holdings. See
## "Reproducibility boundaries" in the top-level README.md.
## Expected input columns: environment, line, block, and the six trait columns.
suppressPackageStartupMessages({library(asreml); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 01_within_environment_blues.R <plot_data.csv> <BLUEs.csv>")
d <- fread(args[1])
traits <- c("DTF", "DTM", "HSW", "PH", "PPP", "YPP")
required <- c("environment", "line", "block", traits)
if (length(miss <- setdiff(required, names(d)))) stop("Missing columns: ", paste(miss, collapse = ", "))
d[, `:=`(environment = factor(environment), line = factor(line), block = factor(block))]

ans <- list()
for (e in levels(d$environment)) {
  de <- droplevels(d[environment == e])
  for (tr in traits) {
    z <- de[is.finite(get(tr)), .(line, block, y = get(tr))]
    fit <- asreml(fixed = y ~ line, random = ~ block,
                  residual = ~ id(units), data = z, trace = FALSE)
    for (k in seq_len(10L)) { if (fit$converge) break; fit <- update(fit) }
    if (!fit$converge) warning("Non-convergence: ", e, " / ", tr)
    pr <- as.data.table(predict(fit, classify = "line", sed = TRUE)$pvals)
    ans[[length(ans) + 1L]] <- pr[, .(environment = e, line,
      Trait = tr, BLUE = predicted.value, SE = standard.error)]
  }
}
fwrite(rbindlist(ans), args[2])
