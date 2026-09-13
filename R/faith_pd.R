# Faith's phylogenetic diversity.
#
# PD is the total branch length of the minimal subtree connecting a sample's
# observed taxa to the root. It is computed by mia::addAlpha(index =
# "faith_diversity"), which the mia documentation states gives values equivalent
# to picante::pd with include.root = TRUE. AmpliPub adds only the guards on the
# tree, which turn a silent wrong answer into a stated refusal.
#
# Cross-checked against `qiime diversity core-metrics-phylogenetic` output on
# the real 487-sample Baxter table; see tests/testthat/test-crosscheck-qiime.R.

#' @keywords internal
ap_faith_pd_mia <- function(counts, tree) {
  ap_assert(inherits(tree, "phylo"), "`tree` must be an `ape::phylo`.")
  ap_assert(!is.null(tree$edge.length),
            paste0("The tree has no branch lengths, so Faith's PD and UniFrac are ",
                   "undefined on it. A cladogram cannot give phylogenetic diversity."))
  ap_assert(all(tree$edge.length >= 0),
            "The tree has {sum(tree$edge.length < 0)} negative branch length{?s}.")
  ap_check_tree_covers(tree, rownames(counts))

  # mia 1.16.1 crashes the R session, rather than raising an error, on a tree
  # whose internal nodes are not numbered in traversal order. ape::root() leaves
  # trees in exactly that state. Writing the tree to Newick and reading it back
  # renumbers the nodes the way ape::read.tree always does, without changing
  # topology or branch lengths (17 significant digits, so nothing is rounded).
  # Trees read from files, including every QIIME tree, are already in this form.
  tree <- ape::read.tree(text = ape::write.tree(tree, digits = 17))

  se <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts), rowTree = tree
  )
  res <- mia::addAlpha(se, index = "faith_diversity")
  stats::setNames(as.numeric(SummarizedExperiment::colData(res)[["faith_diversity"]]),
                  colnames(counts))
}

#' Faith's phylogenetic diversity
#'
#' Total branch length of the minimal subtree connecting each sample's observed
#' features to the root of the phylogeny, computed by
#' `mia::addAlpha(index = "faith_diversity")`.
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
  ap_faith_pd_mia(counts, tree)
}
