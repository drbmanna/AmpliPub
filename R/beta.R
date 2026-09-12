# Beta diversity.
#
# Multiple metrics by default, never one. The four here disagree on purpose:
# each weights abundance and phylogeny differently, and where they disagree is
# informative. Unweighted UniFrac separating groups while weighted UniFrac does
# not says the difference is in which rare taxa are present, not in the
# dominant ones shifting. A paper reporting only Bray-Curtis cannot say that.

#' Beta diversity distances
#'
#' Computes several between-sample distances side by side.
#'
#' | Metric | Weights abundance | Uses the tree |
#' | --- | --- | --- |
#' | `jaccard` | no | no |
#' | `bray_curtis` | yes | no |
#' | `unweighted_unifrac` | no | yes |
#' | `weighted_unifrac` | yes | yes |
#' | `aitchison` | yes, compositionally | no |
#'
#' @section Rarefaction:
#' Bray-Curtis and Jaccard are sensitive to library size, so rarefaction is
#' applied by default, once, with the seed recorded. Aitchison distance works
#' on log-ratios and is scale-invariant, so it does not need it.
#'
#' @section Aitchison distance and zeros:
#' Aitchison distance is Euclidean distance on centred log-ratios, which is the
#' compositionally correct answer and is undefined at zero. The zero
#' replacement is therefore a decision, not a detail: it is done with an
#' explicit pseudocount that is recorded in the result, so a reader can see
#' what was substituted and how much of the table it touched.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param metrics Metrics to compute.
#' @param rarefy Rarefy to a common depth first. Default `TRUE`.
#' @param depth Rarefaction depth. `NULL` uses the smallest library size.
#' @param seed Random seed, recorded with the result.
#' @param pseudocount Value substituted for zeros before the CLR in Aitchison
#'   distance. Default `0.5`.
#'
#' @return An object of class `ap_beta`: a list with `distances` (a named list
#'   of `dist` objects) and the settings used.
#' @export
ap_beta <- function(x,
                    metrics = c("bray_curtis", "jaccard",
                                "unweighted_unifrac", "weighted_unifrac"),
                    rarefy = TRUE,
                    depth = NULL,
                    seed = 1L,
                    pseudocount = 0.5) {

  known <- c("bray_curtis", "jaccard", "unweighted_unifrac", "weighted_unifrac",
             "weighted_normalized_unifrac", "aitchison", "euclidean")
  bad <- setdiff(metrics, known)
  ap_assert(length(bad) == 0L,
            "{cli::qty(length(bad))}Unknown metric{?s}: {paste(bad, collapse = ', ')}. Known: {paste(known, collapse = ', ')}.")

  counts <- SummarizedExperiment::assay(x, "counts")
  ap_check_count_matrix(counts)

  needs_tree <- grepl("unifrac", metrics)
  tree <- tryCatch(TreeSummarizedExperiment::rowTree(x), error = function(e) NULL)
  if (any(needs_tree) && is.null(tree)) {
    ap_warn(paste0(
      "{sum(needs_tree)} metric{?s} need a phylogeny and this object has none. ",
      "Dropping {paste(metrics[needs_tree], collapse = ', ')}."
    ))
    metrics <- metrics[!needs_tree]
  }
  ap_assert(length(metrics) > 0L, "No metrics left to compute.")

  dropped <- character(0)
  if (rarefy) {
    ap_assert(all(counts == round(counts)),
              "Rarefaction needs integer counts. Pass `rarefy = FALSE` for a normalised table.")
    depths <- colSums(counts)
    if (is.null(depth)) {
      depth <- min(depths)
      cli::cli_inform("Rarefying to {format(depth, big.mark = ',')} reads, the smallest library.")
    }
    dropped <- names(depths)[depths < depth]
    if (length(dropped) > 0L) {
      ap_warn("{length(dropped)} sample{?s} below depth {depth} and dropped.")
      counts <- counts[, setdiff(colnames(counts), dropped), drop = FALSE]
    }
    set.seed(seed)
    counts <- ap_rarefy_matrix(counts, depth)
    counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  }

  pd_index <- if (any(grepl("unifrac", metrics))) {
    ap_pd_index(tree, rownames(counts))
  } else NULL

  distances <- list()
  for (m in metrics) {
    distances[[m]] <- switch(
      m,
      bray_curtis = vegan::vegdist(t(counts), method = "bray"),
      jaccard = vegan::vegdist(t(counts), method = "jaccard", binary = TRUE),
      euclidean = stats::dist(t(counts)),
      unweighted_unifrac = ap_unweighted_unifrac(counts, pd_index),
      weighted_unifrac = ap_weighted_unifrac(counts, pd_index, normalized = FALSE),
      weighted_normalized_unifrac = ap_weighted_unifrac(counts, pd_index, normalized = TRUE),
      aitchison = ap_aitchison(counts, pseudocount = pseudocount)
    )
  }

  structure(
    list(distances = distances,
         metrics = metrics,
         rarefied = rarefy,
         depth = if (rarefy) depth else NA_real_,
         seed = if (rarefy) seed else NA_integer_,
         dropped = dropped,
         pseudocount = pseudocount,
         n_samples = ncol(counts),
         metadata = as.data.frame(SummarizedExperiment::colData(x))),
    class = "ap_beta"
  )
}

#' Aitchison distance
#'
#' Euclidean distance between centred log-ratio transformed samples. This is
#' the distance that respects the compositional nature of sequencing data:
#' a feature table records proportions of a fixed read budget, not absolute
#' abundances, and a metric that ignores that reports a rise in one taxon as a
#' fall in every other.
#'
#' The CLR is undefined at zero, and sequencing tables are mostly zeros. The
#' replacement used is recorded in the result rather than applied silently.
#'
#' @param counts Feature-by-sample count matrix.
#' @param pseudocount Value substituted for zeros. Default `0.5`.
#' @return A `stats::dist`.
#' @export
ap_aitchison <- function(counts, pseudocount = 0.5) {
  ap_assert(pseudocount > 0,
            "`pseudocount` must be positive; the CLR is undefined at zero.")
  zero_frac <- mean(counts == 0)
  if (zero_frac > 0.9) {
    ap_warn(paste0(
      "{round(100 * zero_frac, 1)}% of this table is zero, so the Aitchison distance ",
      "is dominated by the pseudocount ({pseudocount}) rather than by the data. ",
      "Filter low-prevalence features first, or prefer a metric defined at zero."
    ))
  }
  m <- counts
  m[m == 0] <- pseudocount
  logm <- log(m)
  clr <- sweep(logm, 2, colMeans(logm), "-")
  stats::dist(t(clr))
}

#' @export
print.ap_beta <- function(x, ...) {
  cli::cli_h1("Beta diversity")
  cli::cli_text("{x$n_samples} samples, {length(x$metrics)} metric{?s}")
  if (x$rarefied) {
    cli::cli_text("Rarefied to {format(x$depth, big.mark = ',')} reads, seed {x$seed}")
  } else {
    cli::cli_alert_warning("Not rarefied. Bray-Curtis and Jaccard track library size.")
  }
  if (length(x$dropped) > 0L) cli::cli_text("{length(x$dropped)} sample{?s} dropped below depth")

  summ <- do.call(rbind, lapply(names(x$distances), function(m) {
    d <- x$distances[[m]]
    data.frame(metric = m, mean = round(mean(d), 4), sd = round(stats::sd(d), 4),
               min = round(min(d), 4), max = round(max(d), 4),
               stringsAsFactors = FALSE)
  }))
  print(summ, row.names = FALSE)
  invisible(x)
}
