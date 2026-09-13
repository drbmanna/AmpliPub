# Unweighted UniFrac.
#
# This is a documented exception to AmpliPub's rule of taking statistics from
# established packages. The only installed package option, rbiom 2.2.1 through
# mia::getDissimilarity, gives a different unweighted answer from scikit-bio
# 0.6.2 (QIIME 2's engine) whenever a pair of samples does not span the root of
# the tree: on four such toy pairs rbiom returned 1/2, 1, 1 and 1/3 where
# scikit-bio and phyloseq return 1/3, 2/3, 1/2 and 1/4. phyloseq agrees with
# scikit-bio on those pairs but differs from QIIME by up to 0.165 on the real
# Baxter table. Weighted UniFrac is taken from rbiom, which agreed everywhere.
#
# The code below is the tip-by-edge incidence formulation: each sample becomes a
# vector over tree edges, and unweighted UniFrac is (union - shared) / union of
# the branch length observed in either sample. It is pinned by the scikit-bio
# reference values in test-beta.R and by QIIME 2 output on 487 real samples in
# test-crosscheck-qiime.R.

#' @keywords internal
ap_pd_index <- function(tree, feature_ids) {
  ap_check_phylo(tree, feature_ids)

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
ap_edge_presence <- function(counts, pd_index) {
  # Edges x samples, 1 where at least one observed tip descends from the edge.
  presence <- Matrix::Matrix((counts > 0) * 1, sparse = TRUE)
  e <- Matrix::t(Matrix::crossprod(presence, pd_index$incidence)) > 0
  as.matrix(e) * 1
}

#' @keywords internal
ap_unweighted_unifrac <- function(counts, pd_index) {
  E <- ap_edge_presence(counts, pd_index)
  b <- pd_index$edge_length

  shared <- crossprod(E, b * E)          # sum_e b_e * E_ei * E_ej
  total <- colSums(b * E)                # branch length observed in each sample
  union <- outer(total, total, "+") - shared

  d <- 1 - shared / union
  d[union == 0] <- 0                     # two empty samples differ by nothing
  diag(d) <- 0
  dimnames(d) <- list(colnames(counts), colnames(counts))
  stats::as.dist(d)
}
