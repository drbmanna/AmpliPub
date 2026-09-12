# Faith's phylogenetic diversity.
#
# PD is the total branch length of the minimal subtree connecting a sample's
# observed taxa to the root. The obvious implementation prunes the tree per
# sample with `ape::keep.tip` and sums the edges. That is correct and far too
# slow here: 490 samples x 100 rarefaction iterations is 49,000 tree prunings.
#
# Instead the tree is reduced once to a tip-by-edge incidence matrix, where
# entry (t, e) is 1 when edge e lies on the path from tip t to the root. An
# edge belongs to a sample's subtree exactly when at least one observed tip sits
# below it, so PD becomes one sparse matrix product per iteration.
#
# Cross-checked against `qiime diversity core-metrics-phylogenetic` output on
# the real 490-sample Baxter table; see tests/testthat/test-crosscheck-qiime.R.

#' @keywords internal
ap_pd_index <- function(tree, feature_ids) {
  ap_assert(inherits(tree, "phylo"), "`tree` must be an `ape::phylo`.")
  ap_check_tree_covers(tree, feature_ids)
  ap_assert(!is.null(tree$edge.length),
            paste0("The tree has no branch lengths, so Faith's PD and UniFrac are ",
                   "undefined on it. A cladogram cannot give phylogenetic diversity."))
  ap_assert(all(tree$edge.length >= 0),
            "The tree has {sum(tree$edge.length < 0)} negative branch length{?s}.")

  n_tip <- length(tree$tip.label)
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  # Edge index by child node: every node except the root is the child of
  # exactly one edge, so this is a complete lookup for walking rootward.
  edge_of_child <- integer(max(tree$edge))
  edge_of_child[child] <- seq_along(child)
  parent_of <- integer(max(tree$edge))
  parent_of[child] <- parent

  i <- integer(0)
  j <- integer(0)
  for (t in seq_len(n_tip)) {
    node <- t
    repeat {
      e <- edge_of_child[node]
      if (e == 0L) break          # reached the root
      i <- c(i, t)
      j <- c(j, e)
      node <- parent_of[node]
    }
  }

  incidence <- Matrix::sparseMatrix(
    i = i, j = j, x = 1,
    dims = c(n_tip, length(child)),
    dimnames = list(tree$tip.label, NULL)
  )

  list(
    incidence = incidence[feature_ids, , drop = FALSE],
    edge_length = tree$edge.length,
    total = sum(tree$edge.length)
  )
}

#' @keywords internal
ap_faith_pd <- function(counts, pd_index) {
  ap_assert(!is.null(pd_index), "Faith's PD was requested without a tree index.")
  ap_assert(
    identical(rownames(counts), rownames(pd_index$incidence)),
    "The tree index was built for a different feature set than this table."
  )
  presence <- Matrix::Matrix((counts > 0) * 1, sparse = TRUE)
  covered <- Matrix::crossprod(presence, pd_index$incidence) > 0
  out <- as.numeric(covered %*% pd_index$edge_length)
  stats::setNames(out, colnames(counts))
}

#' Faith's phylogenetic diversity
#'
#' Total branch length of the minimal subtree connecting each sample's observed
#' features to the root of the phylogeny.
#'
#' Usually called through [ap_alpha()], which rarefies first. PD rises with the
#' number of features observed and therefore with sequencing depth, so raw PD on
#' unequal libraries is not comparable across samples.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()], with a tree.
#' @return A named numeric vector, one value per sample.
#' @export
ap_faith_pd_raw <- function(x) {
  counts <- SummarizedExperiment::assay(x, "counts")
  tree <- tryCatch(TreeSummarizedExperiment::rowTree(x), error = function(e) NULL)
  ap_assert(!is.null(tree),
            "This object has no phylogeny. Pass `tree =` to `ap_import()`.")
  ap_faith_pd(counts, ap_pd_index(tree, rownames(counts)))
}
