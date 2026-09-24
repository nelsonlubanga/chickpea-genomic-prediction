required <- c(
  BGLR = "1.1.4", GAPIT = "4.1.0", asreml = "4.2.0.482",
  data.table = "1.18.4", lme4 = "2.0-6", lmerTest = "3.2-1",
  emmeans = "2.0.4", ggplot2 = "4.0.3", writexl = "2.0.0"
)
status <- lapply(names(required), function(pkg) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  version <- if (installed) as.character(packageVersion(pkg)) else NA_character_
  data.frame(Package = pkg, Required_or_recorded = required[[pkg]],
             Installed = installed, Installed_version = version)
})
status <- do.call(rbind, status)
print(status, row.names = FALSE)
if (!all(status$Installed)) quit(status = 1)
