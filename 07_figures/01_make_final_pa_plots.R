#!/usr/bin/env Rscript

## Final predictive-ability plots used for Figures 2-5 and the combined panel.
## Usage:
##   Rscript 07_figures/01_make_final_pa_plots.R [fold_level_PA.csv] [output_dir]

rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[[1L]] else
  "reference_results/data/PA_fold_level_M1_M4_top500.csv"
output_dir <- if (length(args) >= 2L) args[[2L]] else
  "reference_results/figures"

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pa <- read.csv(input_file, check.names = FALSE)
required <- c("Model", "Panel", "Trait", "CV", "PA")
missing <- setdiff(required, names(pa))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

pa$Trait[pa$Trait == "YPPlnt"] <- "YPP"
pa <- pa %>%
  mutate(
    Model = factor(Model, levels = c("M1", "M2", "M3", "M4")),
    Panel = factor(Panel, levels = c("genome-wide", "haplotype-tagged")),
    Trait = factor(Trait, levels = c("HSW", "PH", "DTF", "DTM", "PPP", "YPP")),
    CV = factor(CV, levels = c("CV2", "CV0", "CV1", "CV00"))
  )

if (anyNA(pa[c("Model", "Panel", "Trait", "CV")])) {
  stop("Unexpected Model, Panel, Trait, or CV label in input data")
}

## One figure per CV scheme (Figures 2-5 layout).
for (cv in levels(pa$CV)) {
  p <- pa %>%
    filter(CV == cv) %>%
    ggplot(aes(x = Trait, y = PA, fill = Model)) +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8)) +
    facet_wrap(~ Panel, nrow = 1) +
    labs(title = cv, x = NULL, y = "Predictive ability", fill = "Model") +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  ggsave(file.path(output_dir, paste0("PA_boxplot_", cv, ".png")),
         p, width = 10, height = 5, dpi = 300)
  ggsave(file.path(output_dir, paste0("PA_boxplot_", cv, ".pdf")),
         p, width = 10, height = 5)
}

## All CV schemes in one trait-by-CV grid.
p_all <- ggplot(pa, aes(x = Model, y = PA, fill = Panel)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  facet_grid(CV ~ Trait) +
  labs(x = NULL, y = "Predictive ability", fill = "Panel") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "PA_boxplot_all.png"),
       p_all, width = 12, height = 9, dpi = 300)
ggsave(file.path(output_dir, "PA_boxplot_all.pdf"),
       p_all, width = 12, height = 9)

cat("Final plots written to", normalizePath(output_dir), "\n")
