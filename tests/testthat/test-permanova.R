# The four cases the dispersion check exists to tell apart. Each is built so the
# right verdict is known before the function runs.

ap_permanova_fixture <- function(counts, meta) {
  x <- ap_import(counts, meta)
  ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
}

test_that("a planted location shift with equal spread is called a location shift", {
  # Both groups are drawn with the same per-sample noise, so their dispersions
  # match by construction; only the underlying proportions differ. Adding a
  # constant to some features instead would shrink one group's relative
  # variability and plant a dispersion difference alongside the shift.
  set.seed(20)
  n_f <- 40L; per <- 25L
  base <- stats::runif(n_f, 5, 50)
  shifted <- base
  shifted[1:8] <- shifted[1:8] * 2.5

  draw <- function(props, n) {
    vapply(seq_len(n), function(i) {
      lambda <- pmax(props * stats::rlnorm(n_f, 0, 0.4), 0.1)
      stats::rmultinom(1, 4000, lambda / sum(lambda))[, 1]
    }, numeric(n_f))
  }
  counts <- cbind(draw(base, per), draw(shifted, per))
  dimnames(counts) <- list(paste0("f", seq_len(n_f)), paste0("s", seq_len(2 * per)))
  meta <- data.frame(g = rep(c("a", "b"), each = per),
                     row.names = colnames(counts), stringsAsFactors = FALSE)

  pn <- ap_permanova(ap_permanova_fixture(counts, meta), "g", permutations = 199L)
  expect_lt(pn$results$p, 0.05)
  expect_gt(pn$dispersion$dispersion_p, 0.05)
  expect_equal(pn$interpretation$verdict, "location shift")
  expect_match(pn$interpretation$interpretation, "genuine shift in centroid")
})

test_that("no difference at all is called no difference", {
  set.seed(21)
  n_f <- 30L; n <- 40L
  counts <- matrix(stats::rpois(n_f * n, 50), nrow = n_f,
                   dimnames = list(paste0("f", seq_len(n_f)), paste0("s", seq_len(n))))
  meta <- data.frame(g = rep(c("a", "b"), each = n / 2),
                     row.names = colnames(counts), stringsAsFactors = FALSE)

  pn <- ap_permanova(ap_permanova_fixture(counts, meta), "g", permutations = 199L)
  expect_gt(pn$results$p, 0.05)
  expect_equal(pn$interpretation$verdict, "no difference")
})

test_that("equal centroids with unequal spread is called dispersion only", {
  # Both groups drawn around the same composition; group b is far noisier. A
  # PERMANOVA alone would report this as null and the real finding would be lost.
  set.seed(22)
  n_f <- 40L; per <- 30L
  base <- stats::runif(n_f, 5, 50)
  make <- function(noise) {
    vapply(seq_len(per), function(i) {
      lambda <- pmax(base * stats::rlnorm(n_f, 0, noise), 0.1)
      stats::rmultinom(1, 4000, lambda / sum(lambda))[, 1]
    }, numeric(n_f))
  }
  counts <- cbind(make(0.15), make(1.2))
  dimnames(counts) <- list(paste0("f", seq_len(n_f)), paste0("s", seq_len(2 * per)))
  meta <- data.frame(g = rep(c("tight", "loose"), each = per),
                     row.names = colnames(counts), stringsAsFactors = FALSE)

  pn <- ap_permanova(ap_permanova_fixture(counts, meta), "g", permutations = 199L)
  expect_lt(pn$dispersion$dispersion_p, 0.05)
  expect_gt(pn$dispersion$max_centroid_ratio, 1.2)
  expect_true(pn$interpretation$verdict %in% c("dispersion only", "confounded by dispersion"))
})

test_that("the dispersion test is always run and always present in the result", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  pn <- ap_permanova(b, "group", permutations = 99L)
  expect_true(all(c("dispersion_F", "dispersion_p", "max_centroid_ratio") %in%
                    names(pn$dispersion)))
  expect_false(is.na(pn$dispersion$dispersion_p[1]))
  expect_equal(nrow(pn$interpretation), nrow(pn$results))
})

test_that("printing shows the verdict alongside the p-value", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  pn <- ap_permanova(b, "group", permutations = 99L)
  # cli writes headings as conditions, the data frame to stdout, so both streams
  # are captured.
  stdout_txt <- paste(utils::capture.output(
    msg_txt <- paste(utils::capture.output(print(pn), type = "message"), collapse = " ")
  ), collapse = " ")
  expect_match(stdout_txt, "verdict")
  expect_match(msg_txt, "What this means", fixed = TRUE)
})

# --- settings that change the answer and are usually left implicit ---

test_that("marginal and sequential sums of squares are both available and by is recorded", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  m <- ap_permanova(b, ~ batch_run + group, permutations = 99L, by = "margin")
  s <- ap_permanova(b, ~ batch_run + group, permutations = 99L, by = "terms")
  expect_equal(m$by, "margin")
  expect_equal(s$by, "terms")
  expect_equal(nrow(m$results), 2L)
})

test_that("the seed makes the permutation p-value reproducible", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  p1 <- ap_permanova(b, "group", permutations = 199L, seed = 42L)
  p2 <- ap_permanova(b, "group", permutations = 199L, seed = 42L)
  expect_equal(p1$results$p, p2$results$p)
  expect_equal(p1$dispersion$dispersion_p, p2$dispersion$dispersion_p)
})

test_that("a continuous term is fitted and marked as having no dispersion test", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  pn <- ap_permanova(b, "age", permutations = 99L)
  expect_equal(pn$results$df, 1L)
  expect_true(is.na(pn$dispersion$dispersion_p))
  expect_match(pn$dispersion$note, "betadisper needs groups")
})

test_that("samples with missing metadata are dropped and the drop is announced", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  cd <- SummarizedExperiment::colData(x)
  cd$group[1:3] <- NA
  SummarizedExperiment::colData(x) <- cd
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  expect_message(pn <- ap_permanova(b, "group", permutations = 99L), "dropped from the")
  expect_equal(pn$results$n, ncol(x) - 3L)
})

test_that("a variable that is not in the metadata is refused", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  expect_error(ap_permanova(b, "nope", permutations = 99L), "not in the sample metadata")
})

test_that("a formula and a character vector of terms give the same result", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  a <- ap_permanova(b, ~ group, permutations = 99L, seed = 3L)
  c2 <- ap_permanova(b, "group", permutations = 99L, seed = 3L)
  expect_equal(a$results$R2, c2$results$R2)
})

# --- plots ---

test_that("the ordination plot builds and carries the PERMANOVA annotation", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  pn <- ap_permanova(b, "group", permutations = 99L)
  ord <- ap_ordinate(b, "bray_curtis")
  p <- ap_plot_ordination(ord, group = "group", permanova = pn)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  expect_match(p$labels$x, "PCo1")
})

test_that("the dispersion plot builds and refuses a continuous term", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  pn <- ap_permanova(b, c("group", "age"), permutations = 99L)
  expect_s3_class(ap_plot_dispersion(pn, "bray_curtis", "group"), "ggplot")
  expect_error(ap_plot_dispersion(pn, "bray_curtis", "age"), "continuous")
})
