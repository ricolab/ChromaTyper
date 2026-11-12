# ChromaTyper: whole-gene chromatin configurations (chromatypes) from ChromHMM states.
# See README.md for the method. Functions only; the command-line entry point is
# chromatyper.R at the repository root.

INVALID_STATES <- c("Type_0", "TYPE_0", "type_0", "0", "", "NA", "N/A", ".")

#' Read an emission matrix: first column the state ID, then one column per mark.
read_emissions <- function(file) {
  e <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  m <- as.matrix(e[-1]); storage.mode(m) <- "numeric"
  rownames(m) <- trimws(e[[1]])
  if (anyNA(m)) stop("Emission matrix has missing or nonnumeric values.")
  if (anyDuplicated(rownames(m))) stop("Duplicated state IDs in the emission matrix.")
  m
}

#' State-to-state distance: sqrt(mean((emission_i - emission_j)^2)).
state_distances <- function(emissions) {
  d <- as.matrix(stats::dist(emissions)) / sqrt(ncol(emissions))
  diag(d) <- 0
  d
}

#' Read the long gene-bin table (CSV, TSV or RDS). bin_weight defaults to 1.
read_gene_bins <- function(file) {
  x <- if (grepl("\\.rds$", file, ignore.case = TRUE)) readRDS(file) else data.table::fread(file)
  x <- data.table::as.data.table(x)
  need <- c("gene_id", "sample_id", "bin_id", "state")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("Gene-bin table lacks column(s): ", paste(miss, collapse = ", "))
  if (!"bin_weight" %in% names(x)) x[, bin_weight := 1]
  x[, c("gene_id", "sample_id", "bin_id", "state") := lapply(.SD, function(v) trimws(as.character(v))),
    .SDcols = need]
  x[state %in% INVALID_STATES, state := NA_character_]
  x
}

#' Sample x sample distance for one gene.
#' states: samples x bins character matrix (NA = missing); w: bin weights.
#' D(A, B) = sqrt(sum(w * d^2) / sum(w)) over bins valid in both samples ("rms"),
#' or sum(w * d) / sum(w) ("mean"). Pairs sharing less than min_valid of the
#' total weight get NA.
gene_distance <- function(states, sd, w = rep(1, ncol(states)),
                          aggregation = c("rms", "mean"), min_valid = 0.8) {
  aggregation <- match.arg(aggregation)
  unk <- setdiff(unique(na.omit(as.vector(states))), rownames(sd))
  if (length(unk)) stop("States absent from the emission matrix: ", paste(unk, collapse = ", "))
  idx <- matrix(match(states, rownames(sd)), nrow(states))
  n <- nrow(states); num <- matrix(0, n, n); den <- matrix(0, n, n)
  for (b in seq_len(ncol(states))) {
    if (w[b] <= 0) next
    ok <- !is.na(idx[, b])
    if (!any(ok)) next
    d <- sd[idx[ok, b], idx[ok, b], drop = FALSE]
    num[ok, ok] <- num[ok, ok] + w[b] * (if (aggregation == "rms") d^2 else d)
    den[ok, ok] <- den[ok, ok] + w[b]
  }
  D <- num / den
  if (aggregation == "rms") D <- sqrt(D)
  D[den / sum(w) < min_valid | den == 0] <- NA
  diag(D) <- 0
  dimnames(D) <- list(rownames(states), rownames(states))
  D
}

#' Chromatypes of one gene from its rows of the gene-bin table.
#' Returns a named character vector (Chromotype_1, ..., or Unassigned),
#' labels numbered in dendrogram order.
chromatype_gene <- function(g, sd, deepSplit = 0, minClusterSize = 5,
                            pamRespectsDendro = FALSE, aggregation = "rms",
                            linkage = "average", min_valid = 0.8,
                            min_bins = 2, min_samples = 3) {
  samples <- unique(g$sample_id); bins <- unique(g$bin_id)
  if (length(samples) < min_samples) stop(length(samples), " samples (< ", min_samples, ")")
  if (length(bins) < min_bins) stop(length(bins), " bins (< ", min_bins, ")")
  w <- g[, .(w = unique(bin_weight)), by = bin_id]
  if (anyDuplicated(w$bin_id)) stop("bin weights differ between samples")
  states <- matrix(NA_character_, length(samples), length(bins), dimnames = list(samples, bins))
  states[cbind(match(g$sample_id, samples), match(g$bin_id, bins))] <- g$state
  D <- gene_distance(states, sd, w$w[match(bins, w$bin_id)], aggregation, min_valid)
  if (anyNA(D)) stop("sample pairs share less than ", min_valid * 100, "% valid bin weight")
  if (all(D == 0)) return(setNames(rep("Chromotype_1", length(samples)), samples))
  hc <- stats::hclust(stats::as.dist(D), method = linkage)
  cl <- dynamicTreeCut::cutreeDynamic(hc, distM = D, method = "hybrid", deepSplit = deepSplit,
                                      minClusterSize = min(minClusterSize, length(samples)),
                                      pamRespectsDendro = pamRespectsDendro, verbose = 0)
  ord <- cl[hc$order]; lev <- unique(ord[ord != 0])
  out <- ifelse(cl == 0, "Unassigned", paste0("Chromotype_", match(cl, lev)))
  setNames(out, samples)
}

#' Chromatypes for every gene in a gene-bin table.
#' Returns list(matrix = sample x gene labels, summary = per-gene counts,
#' failures = genes that could not be processed and why).
chromatyper <- function(bins, emissions, cores = 1, ...) {
  sd <- state_distances(emissions)
  genes <- unique(bins$gene_id); samples <- unique(bins$sample_id)
  by_gene <- split(bins, by = "gene_id", keep.by = TRUE)[genes]
  run <- function(g) tryCatch(chromatype_gene(g, sd, ...), error = function(e) conditionMessage(e))
  res <- if (cores > 1 && .Platform$OS.type != "windows")
    parallel::mclapply(by_gene, run, mc.cores = cores) else lapply(by_gene, run)
  ok <- vapply(res, function(r) !is.null(names(r)), TRUE)
  M <- matrix(NA_character_, length(samples), sum(ok), dimnames = list(samples, genes[ok]))
  for (g in genes[ok]) M[names(res[[g]]), g] <- res[[g]]
  summ <- data.table::data.table(
    gene_id = genes[ok],
    n_samples = vapply(res[ok], length, 1L),
    n_bins = vapply(by_gene[ok], function(g) data.table::uniqueN(g$bin_id), 1L),
    n_chromatypes = vapply(res[ok], function(r) data.table::uniqueN(r[r != "Unassigned"]), 1L),
    n_unassigned = vapply(res[ok], function(r) sum(r == "Unassigned"), 1L))
  fail <- data.table::data.table(gene_id = genes[!ok], error = unlist(res[!ok]))
  list(matrix = M, summary = summ, failures = fail)
}
