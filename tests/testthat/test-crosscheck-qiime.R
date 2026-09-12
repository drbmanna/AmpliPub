# Cross-check against QIIME 2, a system we did not tune.
#
# Self-written tests prove the code does what its author expected. These prove
# it agrees with an independent implementation (scikit-bio, via QIIME 2
# amplicon 2025.7.0) on a real 490-sample dataset. That is worth more than any
# fixture we shaped ourselves, and it is the only check here capable of catching
# an error we built into both the code and its tests.
#
# The reference output is not in this repo, because research data never goes
# into a public package. These tests skip cleanly when it is absent, so CI stays
# green, and they run on the machine where the pipeline was executed.
#
# Point AMPLIPUB_QIIME_DIR at a directory holding a `core_metrics/` subdirectory
# and `tree/rooted_tree.qza` to run them.

ap_qiime_dir <- function() {
  env <- Sys.getenv("AMPLIPUB_QIIME_DIR", unset = "")
  if (nzchar(env) && dir.exists(env)) return(env)
  default <- "//wsl.localhost/Ubuntu-24.04/home/ome/research/baxter2016/q2"
  if (dir.exists(default)) return(default)
  ""
}

skip_without_qiime_output <- function(...) {
  dir <- ap_qiime_dir()
  if (!nzchar(dir)) {
    skip("No QIIME 2 reference output available (set AMPLIPUB_QIIME_DIR).")
  }
  for (f in c(...)) {
    if (!file.exists(file.path(dir, f))) skip(paste("Missing reference file:", f))
  }
  dir
}

test_that("alpha metrics match QIIME 2 on the identical rarefied table", {
  dir <- skip_without_qiime_output(
    "diversity/core_metrics/rarefied_table.qza",
    "diversity/core_metrics/shannon_vector.qza",
    "diversity/core_metrics/observed_features_vector.qza",
    "diversity/core_metrics/evenness_vector.qza",
    "diversity/core_metrics/faith_pd_vector.qza",
    "tree/rooted_tree.qza"
  )
  cm <- file.path(dir, "diversity", "core_metrics")

  # QIIME's own rarefied table, so the comparison tests the metric and not the
  # rarefaction RNG.
  rt <- ap_read_qza(file.path(cm, "rarefied_table.qza"))
  tree <- ap_read_qza(file.path(dir, "tree/rooted_tree.qza"))
  expect_length(unique(colSums(rt)), 1L)

  meta <- data.frame(dummy = rep("a", ncol(rt)), row.names = colnames(rt))
  x <- ap_import(rt, meta, tree = tree)
  a <- ap_alpha(x, metrics = c("q0", "evenness", "faith_pd", "shannon_entropy"),
                rarefy = FALSE)

  pick <- function(metric) {
    v <- a$values[a$values$metric == metric, ]
    stats::setNames(v$value, v$sample_id)
  }
  align <- function(ours, theirs) {
    ids <- intersect(names(ours), names(theirs))
    expect_gt(length(ids), 100L)
    list(ours = unname(ours[ids]), theirs = unname(as.numeric(theirs[ids])))
  }

  # Richness is an integer count. Anything but exact equality is a bug.
  p <- align(pick("q0"), ap_read_qza(file.path(cm, "observed_features_vector.qza")))
  expect_identical(p$ours, p$theirs)

  # QIIME reports Shannon in bits; so does ap_alpha's shannon_entropy.
  p <- align(pick("shannon_entropy"), ap_read_qza(file.path(cm, "shannon_vector.qza")))
  expect_equal(p$ours, p$theirs, tolerance = 1e-10)

  p <- align(pick("evenness"), ap_read_qza(file.path(cm, "evenness_vector.qza")))
  expect_equal(p$ours, p$theirs, tolerance = 1e-10)

  # Faith's PD comes from our own sparse tip-by-edge index, not from a
  # phylogenetics package, so this is the check that it is right. Tolerance is
  # loosened only to absorb the float32 branch lengths in the Newick file.
  p <- align(pick("faith_pd"), ap_read_qza(file.path(cm, "faith_pd_vector.qza")))
  expect_equal(p$ours, p$theirs, tolerance = 1e-6)
})
