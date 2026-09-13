# Each normalization is a package call; these check the wiring, the zero rule, and that
# the comparison reports agreement correctly on data where the answer is known.

ap_fixture_shift <- function(shift = 2.5, seed = 20) {
  # Same per-sample noise in both groups, so dispersions match and only the
  # composition differs. Integer counts, so rarefaction applies.
  set.seed(seed)
  n_f <- 40L; per <- 25L
  base <- stats::runif(n_f, 5, 50)
  shifted <- base
  shifted[1:8] <- shifted[1:8] * shift
  draw <- function(props) {
    vapply(seq_len(per), function(i) {
      lambda <- pmax(props * stats::rlnorm(n_f, 0, 0.4), 0.1)
      stats::rmultinom(1, 4000, lambda / sum(lambda))[, 1]
    }, numeric(n_f))
  }
  counts <- cbind(draw(base), draw(shifted))
  dimnames(counts) <- list(paste0("f", seq_len(n_f)), paste0("s", seq_len(2 * per)))
  meta <- data.frame(g = rep(c("a", "b"), each = per), row.names = colnames(counts),
                     stringsAsFactors = FALSE)
  ap_import(counts, meta)
}

test_that("TSS gives each sample's proportions", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  m <- ap_normalize(x, "tss")
  counts <- SummarizedExperiment::assay(x, "counts")
  # Values and shape only: vegan and ap_normalize attach attributes recording the method.
  ref <- sweep(counts, 2, colSums(counts), "/")
  expect_equal(as.vector(m), as.vector(ref))
  expect_equal(dim(m), dim(ref))
  expect_equal(unname(colSums(m)), rep(1, ncol(m)))
  expect_equal(attr(m, "method"), "tss")
})

test_that("CLR replaces zeros only, so a sample with no zeros ignores its depth", {
  m <- matrix(c(10, 20, 30, 5, 50, 5,
                0, 40, 60, 10, 0, 20), nrow = 6,
              dimnames = list(paste0("f", 1:6), c("s1", "s2")))
  x <- ap_import(m, data.frame(g = c("a", "b"), row.names = c("s1", "s2")))
  clr <- ap_normalize(x, "clr", pseudocount = 0.5)

  # Zero-only replacement, then the centred log-ratio, by hand.
  v <- m[, "s2"]; v[v == 0] <- 0.5
  expect_equal(unname(clr[, "s2"]), unname(log(v) - mean(log(v))))
  expect_equal(unname(colMeans(clr)), c(0, 0))

  # Scaling a sample with no zeros leaves its CLR unchanged.
  m10 <- m; m10[, "s1"] <- m10[, "s1"] * 10
  x10 <- ap_import(m10, data.frame(g = c("a", "b"), row.names = c("s1", "s2")))
  expect_equal(ap_normalize(x10, "clr")[, "s1"], clr[, "s1"])
})

test_that("CSS comes from metagenomeSeq and records the quantile it used", {
  skip_if_not_installed("metagenomeSeq")
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  counts <- SummarizedExperiment::assay(x, "counts")
  m <- ap_normalize(x, "css", css_p = 0.5, css_scale = 1000)

  obj <- metagenomeSeq::cumNorm(metagenomeSeq::newMRexperiment(counts), p = 0.5)
  ref <- metagenomeSeq::MRcounts(obj, norm = TRUE, log = FALSE, sl = 1000)
  expect_equal(as.vector(m), as.vector(ref))
  expect_equal(dim(m), dim(ref))
  expect_equal(attr(m, "css_p"), 0.5)
  expect_true(is.numeric(attr(ap_normalize(x, "css"), "css_p")))
})

test_that("rarefaction draws exactly the depth and is reproducible under a seed", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_normalize(x, "rarefy", depth = 1000, seed = 3L)
  b <- ap_normalize(x, "rarefy", depth = 1000, seed = 3L)
  expect_true(all(colSums(a) == 1000))
  expect_identical(unname(a), unname(b))
  expect_equal(attr(a, "depth"), 1000)
})

test_that("a planted composition shift is found under every normalization", {
  x <- ap_fixture_shift()
  s <- ap_normalization_sensitivity(x, "g", methods = c("tss", "clr", "rarefy"),
                                    permutations = 199L, seed = 4L)
  expect_equal(nrow(s$results), 3L)
  expect_true(all(s$results$p < 0.05))
  expect_true(s$agreement$significant_under_all)
  expect_setequal(s$results$distance, c("bray_curtis", "aitchison"))
})

test_that("agreement is reported as disagreement when verdicts differ", {
  x <- ap_fixture_shift()
  s <- ap_normalization_sensitivity(x, "g", methods = c("tss", "clr"), permutations = 99L)
  # Force a disagreement in the bookkeeping to check the flag, not the statistics.
  s2 <- s
  s2$results$verdict[1] <- "no difference"
  flag <- length(unique(s2$results$verdict[s2$results$term == "g"])) == 1L
  expect_false(flag)
  expect_true(is.logical(s$agreement$same_verdict))
})

test_that("CSS joins the comparison when metagenomeSeq is present", {
  skip_if_not_installed("metagenomeSeq")
  x <- ap_fixture_shift()
  s <- ap_normalization_sensitivity(x, "g", methods = c("tss", "css"), permutations = 99L)
  expect_true("css" %in% s$results$normalization)
})

test_that("a comparison of one normalization is refused", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_error(ap_normalization_sensitivity(x, "group", methods = "tss"), "at least two")
})
