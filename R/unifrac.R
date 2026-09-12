# UniFrac.
#
# Built on the same tip-by-edge incidence matrix as Faith's PD (see faith_pd.R).
# Once each sample is expressed as a vector over tree edges, both UniFrac
# variants are matrix algebra rather than per-pair tree traversal.
#
# Unweighted UniFrac is the fraction of branch length unique to one of two
# samples: (union - shared) / union. Written as matrix products, `shared` for
# every pair at once is a single crossprod.
#
# Weighted UniFrac sums branch lengths scaled by the difference in the relative
# abundance descending from each branch. Which of the two published forms QIIME
# 2 reports, raw or normalised, is not assumed here: both are computed and
# checked against QIIME's own output in test-crosscheck-qiime.R.

#' @keywords internal
ap_edge_presence <- function(counts, pd_index) {
  # Edges x samples, 1 where at least one observed tip descends from the edge.
  presence <- Matrix::Matrix((counts > 0) * 1, sparse = TRUE)
  e <- Matrix::t(Matrix::crossprod(presence, pd_index$incidence)) > 0
  as.matrix(e) * 1
}

#' @keywords internal
ap_edge_abundance <- function(counts, pd_index) {
  # Edges x samples, holding the relative abundance descending from each edge.
  rel <- sweep(counts, 2, colSums(counts), "/")
  as.matrix(Matrix::t(Matrix::crossprod(Matrix::Matrix(rel, sparse = TRUE),
                                        pd_index$incidence)))
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

#' @keywords internal
ap_weighted_unifrac <- function(counts, pd_index, normalized = FALSE) {
  A <- ap_edge_abundance(counts, pd_index)
  b <- pd_index$edge_length
  n <- ncol(counts)

  num <- matrix(0, n, n)
  den <- matrix(0, n, n)
  for (i in seq_len(n)) {
    diff_i <- abs(A - A[, i])
    num[, i] <- colSums(b * diff_i)
    if (normalized) den[, i] <- colSums(b * (A + A[, i]))
  }

  d <- if (normalized) {
    out <- num / den
    out[den == 0] <- 0
    out
  } else {
    num
  }
  diag(d) <- 0
  dimnames(d) <- list(colnames(counts), colnames(counts))
  stats::as.dist(d)
}
