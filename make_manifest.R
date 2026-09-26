#!/usr/bin/env Rscript

root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)))
  commandArgs(trailingOnly = TRUE)[1] else ".", mustWork = TRUE)
files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE,
                    include.dirs = FALSE)
files <- files[!grepl("(^|/)(FILE_MANIFEST.tsv|checksums_md5.tsv)$", files)]
# Local staging copies of the 24 large vector PDFs are excluded. The compact,
# assembled Supplementary Figure S1 and its page manifest are archived instead.
files <- files[!grepl("(^|/)outputs/supplementary_figures/source_partition01_fold1/", files)]
# git's own internal object store is not part of the archive's file inventory.
files <- files[!grepl("(^|/)\\.git(/|$)", files)]
info <- file.info(files)
rel <- substring(files, nchar(root) + 2L)
ext <- tools::file_ext(rel)
manifest <- data.frame(
  File = rel,
  Bytes = info$size,
  Type = ifelse(ext == "", "file", tolower(ext)),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(root, "FILE_MANIFEST.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

md5 <- tools::md5sum(files)
checksums <- data.frame(MD5 = unname(md5), File = rel, stringsAsFactors = FALSE)
write.table(checksums, file.path(root, "checksums_md5.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("Manifested", nrow(manifest), "files\n")
