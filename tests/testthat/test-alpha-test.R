test_that("a planted difference is found with an effect size and an interval", {
  # Group b is shifted by construction, so the test must find it. Anything else
  # means the test is not measuring what it claims.
  set.seed(4)
  n <- 40L
  meta <- data.frame(g = rep(c("a", "b"), each = n / 2),
                     row.names = sprintf("S%02d", seq_len(n)),
                     stringsAsFactors = FALSE)
  vals <- data.frame(
    sample_id = rownames(meta),
    metric = "q1",
    value = c(stats::rnorm(n / 2, 20, 3), stats::rnorm(n / 2, 28, 3)),
    stringsAsFactors = FALSE
  )
  alpha <- structure(list(values = vals, metrics = "q1", rarefied = TRUE,
                          depth = 1000, n_iter = 10L, seed = 1L,
                          dropped = character(0), metadata = meta),
                     class = "ap_alpha")

  res <- ap_alpha_test(alpha, "g", n_boot = 200L)$results
  expect_lt(res$p, 0.001)
  # Roughly 8/3 standard deviations apart, so a large positive effect.
  expect_gt(res$estimate, 1.5)
  expect_lt(res$ci_low, res$estimate)
  expect_gt(res$ci_high, res$estimate)
  expect_false(res$ci_low <= 0 && res$ci_high >= 0)
})

test_that("no difference gives a null result rather than a small p-value", {
  set.seed(6)
  n <- 60L
  meta <- data.frame(g = rep(c("a", "b"), each = n / 2),
                     row.names = sprintf("S%02d", seq_len(n)),
                     stringsAsFactors = FALSE)
  vals <- data.frame(sample_id = rownames(meta), metric = "q1",
                     value = stats::rnorm(n, 20, 3), stringsAsFactors = FALSE)
  alpha <- structure(list(values = vals, metrics = "q1", rarefied = FALSE,
                          depth = NA_real_, n_iter = 1L, seed = NA_integer_,
                          dropped = character(0), metadata = meta),
                     class = "ap_alpha")
  res <- ap_alpha_test(alpha, "g", n_boot = 200L)$results
  expect_gt(res$p, 0.05)
  # A null effect's interval must cover zero.
  expect_lte(res$ci_low, 0)
  expect_gte(res$ci_high, 0)
})

test_that("skewed data routes to the rank-based test, normal data to the parametric one", {
  build <- function(values) {
    n <- length(values)
    meta <- data.frame(g = rep(c("a", "b"), each = n / 2),
                       row.names = sprintf("S%02d", seq_len(n)),
                       stringsAsFactors = FALSE)
    structure(list(values = data.frame(sample_id = rownames(meta), metric = "q1",
                                       value = values, stringsAsFactors = FALSE),
                   metrics = "q1", rarefied = FALSE, depth = NA_real_,
                   n_iter = 1L, seed = NA_integer_, dropped = character(0),
                   metadata = meta),
              class = "ap_alpha")
  }
  set.seed(8)
  normal <- ap_alpha_test(build(stats::rnorm(60, 10, 1)), "g", n_boot = 100L)$results
  expect_equal(normal$test, "Student's t")
  expect_true(normal$assumptions_met)

  skewed <- ap_alpha_test(build(stats::rexp(60, 0.1)), "g", n_boot = 100L)$results
  expect_equal(skewed$test, "Wilcoxon rank-sum")
  expect_false(skewed$assumptions_met)
  expect_equal(skewed$effect, "Cliff's delta")
})

test_that("Hedges' g matches the value worked out by hand", {
  a <- c(10, 12, 14, 16)
  b <- c(20, 22, 24, 26)
  g <- factor(rep(c("a", "b"), each = 4))
  res <- ap_hedges_g(c(a, b), g)
  sp <- sqrt((stats::var(a) + stats::var(b)) / 2)
  d <- (mean(b) - mean(a)) / sp
  expect_equal(res$estimate, d * (1 - 3 / (4 * 8 - 9)))
})

test_that("Cliff's delta is 1 when the groups do not overlap at all", {
  expect_equal(ap_cliffs_point(c(1, 2, 3), c(10, 11, 12)), 1)
  expect_equal(ap_cliffs_point(c(10, 11, 12), c(1, 2, 3)), -1)
  expect_equal(ap_cliffs_point(c(1, 2, 3), c(1, 2, 3)), 0)
})

test_that("three groups route to a three-group test", {
  set.seed(10)
  n <- 60L
  meta <- data.frame(g = rep(c("a", "b", "c"), each = n / 3),
                     row.names = sprintf("S%02d", seq_len(n)),
                     stringsAsFactors = FALSE)
  vals <- data.frame(sample_id = rownames(meta), metric = "q1",
                     value = stats::rnorm(n, 20, 3), stringsAsFactors = FALSE)
  alpha <- structure(list(values = vals, metrics = "q1", rarefied = FALSE,
                          depth = NA_real_, n_iter = 1L, seed = NA_integer_,
                          dropped = character(0), metadata = meta),
                     class = "ap_alpha")
  res <- ap_alpha_test(alpha, "g", n_boot = 100L)$results
  expect_equal(res$k_groups, 3L)
  expect_true(res$test %in% c("One-way ANOVA", "Kruskal-Wallis"))
})

test_that("the Kruskal-Wallis effect is rstatix's eta squared (H), with an unrounded interval", {
  set.seed(14)
  n <- 90L
  meta <- data.frame(g = rep(c("a", "b", "c"), each = n / 3),
                     row.names = sprintf("S%02d", seq_len(n)), stringsAsFactors = FALSE)
  # Skewed values fail the normality check, which routes to the rank-based test.
  values <- stats::rexp(n, 0.1) + rep(c(0, 4, 8), each = n / 3)
  res <- ap_alpha_test(ap_fixture_alpha(values, meta), "g", n_boot = 300L)$results

  expect_equal(res$test, "Kruskal-Wallis")
  expect_equal(res$effect, "eta squared (H)")
  ref <- rstatix::kruskal_effsize(data.frame(value = values, g = factor(meta$g)),
                                  value ~ g)$effsize
  expect_equal(res$estimate, unname(ref))
  # rstatix's own interval is rounded to 0.01. This one comes from boot and
  # must not be.
  expect_false(isTRUE(all.equal(res$ci_low * 100, round(res$ci_low * 100))))
  expect_lte(res$ci_low, res$estimate)
  expect_gte(res$ci_high, res$estimate)
})

test_that("a covariate that explains the difference removes it", {
  # The group difference is entirely driven by the covariate, so adjusting for
  # the covariate must leave nothing behind.
  set.seed(12)
  n <- 80L
  cov <- stats::rnorm(n, 50, 10)
  g <- ifelse(cov > 50, "b", "a")
  meta <- data.frame(g = g, cov = cov, row.names = sprintf("S%02d", seq_len(n)),
                     stringsAsFactors = FALSE)
  vals <- data.frame(sample_id = rownames(meta), metric = "q1",
                     value = 2 * cov + stats::rnorm(n, 0, 2), stringsAsFactors = FALSE)
  alpha <- structure(list(values = vals, metrics = "q1", rarefied = FALSE,
                          depth = NA_real_, n_iter = 1L, seed = NA_integer_,
                          dropped = character(0), metadata = meta),
                     class = "ap_alpha")

  raw <- ap_alpha_test(alpha, "g", n_boot = 100L)$results
  adj <- ap_alpha_test(alpha, "g", covariates = "cov")$results
  # Unadjusted: a huge apparent group effect, entirely borrowed from cov.
  expect_lt(raw$p, 0.001)
  expect_gt(raw$estimate, 2)
  # Adjusted: gone. A little partial eta squared survives because g is a hard
  # threshold on cov, and a linear term cannot absorb a step exactly.
  expect_gt(adj$p, 0.05)
  expect_lt(adj$estimate, 0.05)
})

test_that("a subject variable routes to a mixed model", {
  skip_if_not_installed("lme4")
  set.seed(14)
  n_subj <- 20L
  subj <- rep(sprintf("s%02d", seq_len(n_subj)), each = 2)
  g <- rep(c("pre", "post"), times = n_subj)
  offset <- rep(stats::rnorm(n_subj, 0, 5), each = 2)
  meta <- data.frame(g = g, subj = subj,
                     row.names = sprintf("S%02d", seq_along(g)),
                     stringsAsFactors = FALSE)
  vals <- data.frame(sample_id = rownames(meta), metric = "q1",
                     value = 20 + offset + ifelse(g == "post", 3, 0) + stats::rnorm(length(g), 0, 1),
                     stringsAsFactors = FALSE)
  alpha <- structure(list(values = vals, metrics = "q1", rarefied = FALSE,
                          depth = NA_real_, n_iter = 1L, seed = NA_integer_,
                          dropped = character(0), metadata = meta),
                     class = "ap_alpha")

  res <- ap_alpha_test(alpha, "g", subject = "subj")$results
  expect_match(res$test, "mixed model")
  # The paired signal is strong once the subject offsets are accounted for.
  expect_lt(res$p, 0.001)
})

test_that("naming a variable that is not there is refused", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = "q0", n_iter = 2L)
  expect_error(ap_alpha_test(a, "nope"), "not in the sample metadata")
})

test_that("p-values across metrics are adjusted for multiplicity", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = c("q0", "q1", "q2"), n_iter = 3L)
  res <- ap_alpha_test(a, "group", n_boot = 100L)$results
  expect_true(all(res$p_adj >= res$p))
  expect_equal(res$p_adj, stats::p.adjust(res$p, method = "BH"))
})

# --- plots ---

test_that("the alpha plot builds with the expected layers", {
  x <- ap_fixture_object(tree = TRUE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = c("q0", "faith_pd"), n_iter = 3L)
  p <- ap_plot_alpha(a, "group", annotate = FALSE)
  expect_s3_class(p, "ggplot")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomBoxplot"), logical(1))))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("the annotated plot carries the test result as a text layer", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = "q0", n_iter = 3L)
  p <- ap_plot_alpha(a, "group")
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1))))
})

test_that("the violin variant builds", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = "q1", n_iter = 3L)
  p <- ap_plot_alpha(a, "group", type = "violin", annotate = FALSE)
  expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomViolin"), logical(1))))
})

test_that("plotting by a variable that is not there is refused", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = "q0", n_iter = 2L)
  expect_error(ap_plot_alpha(a, "nope"), "not in the sample metadata")
})

test_that("the rarefaction curve builds and rises with depth", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  p <- ap_plot_rarefaction(x, depths = c(200, 800, 2000), n_iter = 2L)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  d <- built$data[[1]]
  means <- tapply(d$y, d$x, mean)
  expect_true(all(diff(means) > 0))
})
