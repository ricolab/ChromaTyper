#!/usr/bin/env Rscript
# Usage:
#   Rscript chromatyper.R --emissions emissions.csv --bins gene_bins.csv --out results/prefix
#     [--deepSplit 0] [--minClusterSize 5] [--pamRespectsDendro FALSE]
#     [--aggregation rms] [--minValid 0.8] [--cores 1]
# Writes <out>_chromatypes.csv (sample x gene), <out>_summary.csv, <out>_failures.csv.

suppressPackageStartupMessages(library(data.table))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "R", "chromatyper.R"))

a <- commandArgs(TRUE)
arg <- function(name, default = NULL) {
  i <- which(a == paste0("--", name))
  if (length(i)) a[i + 1] else if (is.null(default)) stop("--", name, " is required") else default
}
emissions <- read_emissions(arg("emissions"))
bins <- read_gene_bins(arg("bins"))
out <- arg("out")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

res <- chromatyper(bins, emissions,
  cores = as.integer(arg("cores", "1")),
  deepSplit = as.integer(arg("deepSplit", "0")),
  minClusterSize = as.integer(arg("minClusterSize", "5")),
  pamRespectsDendro = as.logical(arg("pamRespectsDendro", "FALSE")),
  aggregation = arg("aggregation", "rms"),
  min_valid = as.numeric(arg("minValid", "0.8")))

fwrite(as.data.table(res$matrix, keep.rownames = "sample_id"), paste0(out, "_chromatypes.csv"))
fwrite(res$summary, paste0(out, "_summary.csv"))
fwrite(res$failures, paste0(out, "_failures.csv"))
cat(sprintf("%d genes done, %d failed; output: %s_*.csv\n", nrow(res$summary), nrow(res$failures), out))
