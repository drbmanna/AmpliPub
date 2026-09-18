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
#' Bray-Curtis and Jaccard come from `vegan::vegdist`, weighted UniFrac from
#' `mia::getDissimilarity(method = "unifrac")` (which calls `rbiom::unifrac`),
#' and Aitchison from `vegan::vegdist(method = "aitchison")`. Weighted UniFrac is
#' the raw form, the one QIIME 2 reports; rbiom offers no normalised form.
#'
#' Unweighted UniFrac is computed by AmpliPub itself, the one exception. rbiom
#' gives a different unweighted value from scikit-bio (QIIME 2's engine) for any
#' pair of samples that does not span the root of the tree, and phyloseq differs
#' from QIIME 2 on real data. AmpliPub's version is tested against scikit-bio
#' reference values and against QIIME 2 output.
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
             "aitchison", "euclidean")
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

  distances <- list()
  for (m in metrics) {
    distances[[m]] <- switch(
      m,
      bray_curtis = vegan::vegdist(t(counts), method = "bray"),
      jaccard = vegan::vegdist(t(counts), method = "jaccard", binary = TRUE),
      euclidean = stats::dist(t(counts)),
      unweighted_unifrac = ap_unifrac(counts, tree, weighted = FALSE),
      weighted_unifrac = ap_unifrac(counts, tree, weighted = TRUE),
      aitchison = ap_aitchison(counts, pseudocount = pseudocount)
    )
  }

  structure(
    list(distances = distances,
         metrics = metrics,
         rarefied = rarefy,
         depth = if (rarefy) depth else NA_real_,
         common_depth = ap_common_depth(colSums(counts)),
         seed = if (rarefy) seed else NA_integer_,
         dropped = dropped,
         pseudocount = pseudocount,
         n_samples = ncol(counts),
         metadata = as.data.frame(SummarizedExperiment::colData(x))),
    class = "ap_beta"
  )
}

# UniFrac. Weighted comes from mia::getDissimilarity(method = "unifrac"), which
# calls rbiom::unifrac and computes the raw form QIIME 2 reports. Unweighted comes
# from AmpliPub's own code in unifrac.R, because rbiom's unweighted value departs
# from scikit-bio for pairs that do not span the root (see that file). Both keep
# the sample order of the table.
#' @keywords internal
ap_unifrac <- function(counts, tree, weighted) {
  ap_check_phylo(tree, rownames(counts))
  if (!weighted) {
    return(ap_unweighted_unifrac(counts, ap_pd_index(tree, rownames(counts))))
  }
  se <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts), rowTree = tree
  )
  d <- mia::getDissimilarity(se, method = "unifrac", weighted = weighted, tree = tree)
  ids <- colnames(counts)
  stats::as.dist(as.matrix(d)[ids, ids])
}

# The guards every phylogenetic metric needs. A tree without branch lengths, or
# with negative ones, gives a number that looks like diversity and is not.
#' @keywords internal
ap_check_phylo <- function(tree, feature_ids) {
  ap_assert(inherits(tree, "phylo"), "`tree` must be an `ape::phylo`.")
  ap_assert(!is.null(tree$edge.length),
            paste0("The tree has no branch lengths, so Faith's PD and UniFrac are ",
                   "undefined on it. A cladogram cannot give phylogenetic diversity."))
  ap_assert(all(tree$edge.length >= 0),
            "The tree has {sum(tree$edge.length < 0)} negative branch length{?s}.")
  ap_check_tree_covers(tree, feature_ids)
  invisible(TRUE)
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
#' Zeros, and only zeros, are replaced by `pseudocount`, and the distance is then
#' computed by `vegan::vegdist(method = "aitchison")` on the now-positive table.
#' vegan's own `pseudocount` argument is not used because it adds the value to
#' every count, which makes the distance change when a sample's depth is scaled.
#' Scale invariance is the property Aitchison distance exists to have.
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
  vegan::vegdist(t(m), method = "aitchison")
}

#' @export
print.ap_beta <- function(x, ...) {
  cli::cli_h1("Beta diversity")
  cli::cli_text("{x$n_samples} samples, {length(x$metrics)} metric{?s}")
  if (x$rarefied) {
    cli::cli_text("Rarefied to {format(x$depth, big.mark = ',')} reads, seed {x$seed}")
  } else if (!is.na(x$common_depth %||% NA_real_)) {
    cli::cli_text("All samples at a common depth of {format(x$common_depth, big.mark = ',')} reads")
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
