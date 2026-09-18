# Every output that a reader has to interpret must carry text explaining it, and
# that text must exist for every combination of inputs, not just the ones the
# Baxter data happened to produce. These tests walk each branch and then assert,
# over a grid, that no combination falls through to an empty explanation.

# --- ordination diagnostics ---------------------------------------------------

ap_fake_pcoa <- function(neg, ax12) {
  structure(
    list(method = "pcoa", metric = "bray_curtis", k = 2L, seed = 1L,
         n_samples = 50L, stress = NA_real_,
         negative_eigenvalue_fraction = neg,
         prop_explained = c(ax12 * 0.6, ax12 * 0.4, rep(0.01, 5))),
    class = "ap_ordination"
  )
}

ap_fake_nmds <- function(stress, converged = TRUE) {
  structure(
    list(method = "nmds", metric = "bray_curtis", k = 2L, seed = 1L,
         n_samples = 50L, stress = stress, converged = converged,
         negative_eigenvalue_fraction = NA_real_, prop_explained = NULL),
    class = "ap_ordination"
  )
}

test_that("PCoA distortion and low variance are each reported, and together", {
  expect_equal(ap_ordination_interpret(ap_fake_pcoa(0.01, 0.55))$verdict, "readable")
  expect_equal(ap_ordination_interpret(ap_fake_pcoa(0.20, 0.55))$verdict, "distorted projection")
  expect_equal(ap_ordination_interpret(ap_fake_pcoa(0.01, 0.10))$verdict,
               "low variance on the plotted axes")
  expect_equal(ap_ordination_interpret(ap_fake_pcoa(0.20, 0.10))$verdict,
               "distorted and low-variance")
})

test_that("the PCoA text names the actual numbers, not a generic caution", {
  it <- ap_ordination_interpret(ap_fake_pcoa(0.20, 0.10))$interpretation
  expect_true(grepl("20.0%", it, fixed = TRUE))
  expect_true(grepl("10.0%", it, fixed = TRUE))
})

test_that("NMDS stress bands are distinguished at the conventional cut points", {
  expect_equal(ap_ordination_interpret(ap_fake_nmds(0.05))$verdict, "good")
  expect_equal(ap_ordination_interpret(ap_fake_nmds(0.15))$verdict, "usable")
  expect_equal(ap_ordination_interpret(ap_fake_nmds(0.25))$verdict,
               "too high to read as a map")
})

test_that("a stress of exactly 0.2 is not called usable", {
  expect_equal(ap_ordination_interpret(ap_fake_nmds(0.2))$verdict,
               "too high to read as a map")
})

test_that("failure to converge is added to the verdict, not substituted for it", {
  d <- ap_ordination_interpret(ap_fake_nmds(0.05, converged = FALSE))
  expect_true(grepl("did not converge", d$verdict, fixed = TRUE))
  expect_true(grepl("good", d$verdict, fixed = TRUE))
  expect_true(grepl("trymax", d$interpretation, fixed = TRUE))
})

test_that("a missing stress value is reported rather than treated as good", {
  expect_equal(ap_ordination_interpret(ap_fake_nmds(NA_real_))$verdict, "not available")
})

test_that("diagnostics accept a list and return one row per ordination", {
  d <- ap_ordination_diagnostics(list(ap_fake_pcoa(0.01, 0.5), ap_fake_nmds(0.15)))
  expect_equal(nrow(d), 2L)
  expect_false(any(is.na(d$interpretation)))
})

test_that("no combination of PCoA or NMDS inputs falls through without text", {
  grid <- expand.grid(neg = c(0, 0.04, 0.05, 0.06, 0.5),
                      ax12 = c(0.01, 0.19, 0.20, 0.21, 0.99))
  for (i in seq_len(nrow(grid))) {
    d <- ap_ordination_interpret(ap_fake_pcoa(grid$neg[i], grid$ax12[i]))
    expect_true(!is.na(d$verdict) && nzchar(d$interpretation))
  }
  for (s in c(0, 0.099, 0.1, 0.199, 0.2, 0.9)) {
    for (cv in c(TRUE, FALSE)) {
      d <- ap_ordination_interpret(ap_fake_nmds(s, cv))
      expect_true(!is.na(d$verdict) && nzchar(d$interpretation))
    }
  }
})

# --- alpha diversity ----------------------------------------------------------

ap_fake_alpha_test <- function(p_adj, assumptions_met, ci = TRUE) {
  res <- data.frame(
    metric = "q1", test = "t-test", n = 40L, k_groups = 2L, statistic = 2,
    p = p_adj, p_adj = p_adj, effect = "Hedges g", estimate = 0.42,
    ci_low = if (ci) 0.1 else NA_real_, ci_high = if (ci) 0.75 else NA_real_,
    normality_p = NA_real_, variance_p = NA_real_,
    assumptions_met = assumptions_met, stringsAsFactors = FALSE
  )
  list(results = res)
}

test_that("a significant alpha result says which test assumption state produced it", {
  expect_equal(ap_alpha_test_interpret(ap_fake_alpha_test(0.01, TRUE))$verdict, "difference")
  expect_equal(ap_alpha_test_interpret(ap_fake_alpha_test(0.01, FALSE))$verdict,
               "difference, rank-based test")
  expect_equal(ap_alpha_test_interpret(ap_fake_alpha_test(0.01, NA))$verdict, "difference")
})

test_that("a rank-based result is not described as a statement about means", {
  it <- ap_alpha_test_interpret(ap_fake_alpha_test(0.01, FALSE))$interpretation
  expect_true(grepl("about ranks, not about means", it, fixed = TRUE))
})

test_that("a null alpha result is reported as absence of evidence, with the interval", {
  it <- ap_alpha_test_interpret(ap_fake_alpha_test(0.4, TRUE))$interpretation
  expect_true(grepl("absence of evidence", it, fixed = TRUE))
  expect_true(grepl("[0.100, 0.750]", it, fixed = TRUE))
})

test_that("a null alpha result without an interval says so instead of implying absence", {
  it <- ap_alpha_test_interpret(ap_fake_alpha_test(0.4, TRUE, ci = FALSE))$interpretation
  expect_true(grepl("not that one is absent", it, fixed = TRUE))
})

test_that("a metric with no p-value is reported as not testable", {
  expect_equal(ap_alpha_test_interpret(ap_fake_alpha_test(NA_real_, NA))$verdict,
               "not testable")
})

test_that("no combination of alpha inputs falls through without text", {
  for (p in c(NA_real_, 0, 0.049, 0.05, 0.051, 1)) {
    for (a in list(TRUE, FALSE, NA)) {
      for (ci in c(TRUE, FALSE)) {
        it <- ap_alpha_test_interpret(ap_fake_alpha_test(p, a[[1]], ci))
        expect_true(!is.na(it$verdict) && nzchar(it$interpretation))
      }
    }
  }
})

# --- normalization sensitivity ------------------------------------------------

ap_fake_norm <- function(p, verdict) {
  data.frame(normalization = c("tss", "css", "clr", "rarefy"), p = p,
             verdict = verdict, stringsAsFactors = FALSE)
}

test_that("a conclusion stable across normalizations is named as the strongest form", {
  it <- ap_normalization_interpret(
    ap_fake_norm(rep(0.001, 4), rep("location shift", 4)), "dx")
  expect_true(grepl("does not depend on the normalization", it, fixed = TRUE))
})

test_that("a stable null is reported as stable, not as a failure", {
  it <- ap_normalization_interpret(ap_fake_norm(rep(0.6, 4), rep("no difference", 4)), "dx")
  expect_true(grepl("null result is stable", it, fixed = TRUE))
})

test_that("agreement on significance but not on meaning is reported as such", {
  it <- ap_normalization_interpret(
    ap_fake_norm(rep(0.001, 4),
                 c("location shift", "confounded by dispersion",
                   "location shift", "location shift")), "dx")
  expect_true(grepl("do not agree on what", it, fixed = TRUE))
})

test_that("a result that depends on the normalization names which ones disagreed", {
  it <- ap_normalization_interpret(
    ap_fake_norm(c(0.001, 0.001, 0.6, 0.6), rep("location shift", 4)), "dx")
  expect_true(grepl("significant under tss, css", it, fixed = TRUE))
  expect_true(grepl("not under clr, rarefy", it, fixed = TRUE))
})

test_that("no combination of normalization outcomes falls through without text", {
  ps <- list(rep(0.001, 4), rep(0.6, 4), c(0.001, 0.6, 0.6, 0.6), c(0.001, 0.001, 0.001, 0.6))
  vs <- list(rep("location shift", 4), c("a", "b", "c", "d"))
  for (p in ps) {
    for (v in vs) {
      it <- ap_normalization_interpret(ap_fake_norm(p, v), "dx")
      expect_true(nzchar(it) && !is.na(it))
    }
  }
})

# --- capturing what a print method emits --------------------------------------
#
# The report is built from these captures. When this helper is wrong the report
# silently loses text rather than failing, which is how every report before
# 2026-09-18 came to show raw list internals.

test_that("cli-only output is captured", {
  meta <- data.frame(dx = rep(c("case", "control"), each = 6),
                     row.names = sprintf("S%02d", 1:12), stringsAsFactors = FALSE)
  scan <- ap_scan_metadata(meta, group = "dx")
  out <- paste(ap_capture_print(scan), collapse = "\n")
  expect_true(grepl("Metadata assessment", out, fixed = TRUE))
  expect_true(grepl("Grouping variable", out, fixed = TRUE))
})

test_that("base print output is captured", {
  out <- paste(ap_capture_print(data.frame(a = 1:2, b = c("x", "y"))), collapse = "\n")
  expect_true(grepl("a b", out, fixed = TRUE))
  expect_true(grepl("1 x", out, fixed = TRUE))
})

test_that("cli and base output keep the order the method emitted them in", {
  obj <- structure(list(), class = "ap_order_probe")
  registerS3method("print", "ap_order_probe", function(x, ...) {
    cli::cli_text("FIRST_CLI")
    print(data.frame(MIDDLE_BASE = 1L), row.names = FALSE)
    cli::cli_text("LAST_CLI")
    invisible(x)
  }, envir = globalenv())
  on.exit(rm("print.ap_order_probe", envir = globalenv()), add = TRUE)

  out <- paste(ap_capture_print(obj), collapse = "\n")
  expect_true(all(c(regexpr("FIRST_CLI", out, fixed = TRUE),
                    regexpr("MIDDLE_BASE", out, fixed = TRUE),
                    regexpr("LAST_CLI", out, fixed = TRUE)) > 0))
  expect_lt(regexpr("FIRST_CLI", out, fixed = TRUE),
            regexpr("MIDDLE_BASE", out, fixed = TRUE))
  expect_lt(regexpr("MIDDLE_BASE", out, fixed = TRUE),
            regexpr("LAST_CLI", out, fixed = TRUE))
})

test_that("the sink is released when the print method errors", {
  obj <- structure(list(), class = "ap_boom_probe")
  registerS3method("print", "ap_boom_probe", function(x, ...) stop("boom"),
                   envir = globalenv())
  on.exit(rm("print.ap_boom_probe", envir = globalenv()), add = TRUE)

  before <- sink.number()
  expect_error(ap_capture_print(obj), "boom")
  expect_equal(sink.number(), before)
})

test_that("the direction-disagreement warning fires when methods disagree on sign", {
  # Baxter produced zero disagreements, so this branch is never exercised by the
  # real data. It is exercised here instead.
  feats <- data.frame(
    feature = c("f1", "f2"), n_methods = c(2L, 1L),
    methods = c("a,b", "a"), direction_agrees = c(FALSE, TRUE),
    consensus = c(FALSE, FALSE), stringsAsFactors = FALSE
  )
  obj <- structure(
    list(features = feats, methods = c("a", "b"), alpha = 0.05, min_methods = 2L,
         contrast = "case_vs_control", n_by_method = c(a = 2L, b = 1L),
         pairwise = matrix(1, 2, 2, dimnames = list(c("a", "b"), c("a", "b"))),
         consensus = character(0)),
    class = "ap_da_concordance"
  )
  out <- paste(ap_capture_print(obj), collapse = " ")
  expect_true(grepl("opposite signs", out, fixed = TRUE))
  expect_true(grepl("1 feature", out, fixed = TRUE))
})
