# Differential abundance is tested against a planted truth, with the success
# criterion written down before the result is looked at. Post hoc, any list of
# significant features looks like a finding.

ap_da_planted <- function(n_features = 60L, per = 20L, n_planted = 6L,
                          fold = 5, depth = 20000, seed = 101) {
  set.seed(seed)
  ids <- paste0("f", sprintf("%02d", seq_len(n_features)))
  planted <- ids[seq_len(n_planted)]
  base <- stats::runif(n_features, 2, 200)
  names(base) <- ids

  draw <- function(props, n) {
    vapply(seq_len(n), function(i) {
      lambda <- pmax(props * stats::rlnorm(n_features, 0, 0.5), 0.05)
      stats::rmultinom(1, depth, lambda / sum(lambda))[, 1]
    }, numeric(n_features))
  }
  case <- base
  case[planted] <- case[planted] * fold

  counts <- cbind(draw(base, per), draw(case, per))
  dimnames(counts) <- list(ids, sprintf("s%02d", seq_len(2 * per)))
  meta <- data.frame(
    dx = rep(c("control", "case"), each = per),
    age = round(stats::runif(2 * per, 30, 70)),
    row.names = colnames(counts), stringsAsFactors = FALSE
  )
  list(x = ap_import(counts, meta), planted = planted, ids = ids)
}

test_that("the bias-corrected methods recover the planted features and stay specific", {
  skip_on_cran()
  skip_if_no_ancombc2()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = c("ancombc2", "aldex2", "linda"),
                               reference = "control", seed = 1L, mc_samples = 64L))
  r <- da$results
  n_null <- length(fx$ids) - length(fx$planted)

  for (m in unique(r$method)) {
    called <- r$feature[r$method == m & r$significant]
    recovered <- intersect(called, fx$planted)
    false_pos <- setdiff(called, fx$planted)
    # Criterion fixed in advance: most of the truth, few of the nulls.
    expect_gte(length(recovered), 4L)
    expect_lt(length(false_pos), 0.15 * n_null)
    # Planted features were enriched in "case", so every recovered effect is
    # positive against the "control" reference.
    expect_true(all(r$effect[r$method == m & r$feature %in% recovered] > 0))
  }
})

test_that("MaAsLin2 on TSS shows compositional bias, and the bias is one-directional", {
  # Not a bug in the adapter, and worth a test because it is the reason the
  # package runs four methods. Boosting some features depresses every other
  # feature's relative abundance; a method without bias correction reads that
  # depression as real depletion. The signature is that the extra calls are all
  # negative.
  skip_on_cran()
  skip_if_not_installed("Maaslin2")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = "maaslin2",
                               reference = "control", seed = 1L))
  r <- da$results
  false_pos <- r[r$significant & !(r$feature %in% fx$planted), ]
  skip_if(nrow(false_pos) == 0L, "No compositional false positives at this seed.")
  expect_true(all(false_pos$effect < 0))
  expect_true(all(r$effect[r$significant & r$feature %in% fx$planted] > 0))
})

test_that("the same prevalence filter reaches every method, so denominators match", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = c("aldex2", "linda"),
                               reference = "control", prv_cut = 0.5, seed = 1L,
                               mc_samples = 32L))
  per_method <- table(da$results$method)
  expect_equal(length(unique(as.integer(per_method))), 1L)
  expect_equal(unique(as.integer(per_method)), da$n_features)
})

test_that("every method reports the same contrast label so they can be compared", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = c("aldex2", "linda"),
                               reference = "control", seed = 1L, mc_samples = 32L))
  expect_equal(unique(da$results$contrast), "case_vs_control")
})

test_that("effect scales are carried per method rather than silently unified", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = c("aldex2", "linda"),
                               reference = "control", seed = 1L, mc_samples = 32L))
  scales <- unique(da$results[, c("method", "effect_scale")])
  expect_equal(nrow(scales), 2L)
  expect_equal(length(unique(scales$effect_scale)), 2L)
})

test_that("reversing the reference level reverses the sign of every effect", {
  skip_on_cran()
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  a <- suppressMessages(ap_da(fx$x, "dx", method = "linda", reference = "control", seed = 1L))
  b <- suppressMessages(ap_da(fx$x, "dx", method = "linda", reference = "case", seed = 1L))
  ea <- a$results$effect[match(fx$ids, a$results$feature)]
  eb <- b$results$effect[match(fx$ids, b$results$feature)]
  ok <- !is.na(ea) & !is.na(eb)
  expect_equal(ea[ok], -eb[ok], tolerance = 1e-6)
})

# --- guards ---

test_that("an already-normalised table is refused", {
  fx <- ap_da_planted(n_features = 20L, per = 6L)
  counts <- SummarizedExperiment::assay(fx$x, "counts")
  rel <- sweep(counts, 2, colSums(counts), "/")
  meta <- as.data.frame(SummarizedExperiment::colData(fx$x))
  y <- suppressWarnings(ap_import(rel, meta, expect_counts = FALSE))
  expect_error(suppressMessages(ap_da(y, "dx", method = "linda")), "already been normalised")
})

test_that("a continuous grouping variable is refused with the reason", {
  fx <- ap_da_planted(n_features = 20L, per = 6L)
  expect_error(suppressMessages(ap_da(fx$x, "age", method = "linda")), "continuous")
})

test_that("an unknown method is refused by name", {
  fx <- ap_da_planted(n_features = 20L, per = 6L)
  expect_error(ap_da(fx$x, "dx", method = "deseq2"), "deseq2")
})

test_that("a prevalence filter that leaves nothing is refused", {
  fx <- ap_da_planted(n_features = 20L, per = 6L)
  expect_error(suppressMessages(ap_da(fx$x, "dx", method = "linda", prv_cut = 1.01)),
               "Lower `prv_cut`")
})

test_that("ALDEx2 refuses covariates rather than quietly ignoring them", {
  skip_if_not_installed("ALDEx2")
  fx <- ap_da_planted(n_features = 20L, per = 8L)
  # It was the only method asked for, so the run cannot continue. The abort has
  # to carry ALDEx2's own reason, not just say that everything failed.
  expect_error(
    suppressWarnings(suppressMessages(
      ap_da(fx$x, "dx", method = "aldex2", covariates = "age",
            seed = 1L, mc_samples = 16L)
    )),
    "takes no covariates"
  )
})

test_that("one failing method does not stop the others, and is reported as skipped", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted(n_features = 30L, per = 10L)
  expect_warning(
    da <- suppressMessages(ap_da(fx$x, "dx", method = c("aldex2", "linda"),
                                 covariates = "age", seed = 1L, mc_samples = 16L)),
    "takes no covariates"
  )
  expect_equal(da$skipped, "aldex2")
  expect_equal(da$methods, "linda")
  expect_gt(nrow(da$results), 0L)
})

test_that("a bad reference level is refused, listing the real ones", {
  fx <- ap_da_planted(n_features = 20L, per = 6L)
  expect_error(suppressMessages(ap_da(fx$x, "dx", method = "linda", reference = "healthy")),
               "Levels: case, control|not a level")
})

# --- concordance ---

test_that("the consensus set is exactly the planted features", {
  skip_on_cran()
  skip_if_no_ancombc2()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = c("ancombc2", "aldex2", "linda"),
                               reference = "control", seed = 1L, mc_samples = 64L))
  cc <- ap_da_concordance(da)
  expect_setequal(cc$consensus, fx$planted)
  expect_true(all(cc$features$direction_agrees))
  expect_equal(dim(cc$pairwise), c(3L, 3L))
  expect_equal(unname(diag(cc$pairwise)), rep(1, 3))
})

test_that("lowering min_methods widens the consensus set", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = c("aldex2", "linda"),
                               reference = "control", seed = 1L, mc_samples = 64L))
  strict <- ap_da_concordance(da, min_methods = 2L)
  loose <- ap_da_concordance(da, min_methods = 1L)
  expect_gte(length(loose$consensus), length(strict$consensus))
})

test_that("a feature no method called is not in the consensus set", {
  skip_on_cran()
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = "linda", reference = "control", seed = 1L))
  cc <- ap_da_concordance(da)
  never <- cc$features$feature[cc$features$n_methods == 0L]
  expect_length(intersect(never, cc$consensus), 0L)
})

# --- plots ---

test_that("the DA plots build", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted()
  da <- suppressMessages(ap_da(fx$x, "dx", method = c("aldex2", "linda"),
                               reference = "control", seed = 1L, mc_samples = 32L))
  cc <- ap_da_concordance(da)

  expect_s3_class(ap_plot_volcano(da), "ggplot")
  expect_s3_class(ap_plot_concordance(cc), "ggplot")
  if (length(cc$consensus) > 0L) {
    expect_s3_class(ap_plot_da_effects(cc), "ggplot")
  }
})

test_that("plotting an empty consensus set says how to widen it", {
  skip_on_cran()
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted(fold = 1, n_planted = 1L)
  da <- suppressMessages(ap_da(fx$x, "dx", method = "linda", reference = "control", seed = 1L))
  cc <- ap_da_concordance(da)
  skip_if(length(cc$consensus) > 0L, "This seed produced calls on a null dataset.")
  expect_error(ap_plot_da_effects(cc), "Lower `min_methods`")
})

test_that("feature IDs survive every backend unchanged, digits and all", {
  # Maaslin2 runs feature names through make.names(), which prepends "X" to any
  # name starting with a digit. ASV hashes start with a digit roughly 40% of the
  # time. The result still looks complete; it just no longer joins to the other
  # methods. This test uses IDs that trip it.
  skip_on_cran()
  skip_if_not_installed("Maaslin2")
  skip_if_not_installed("MicrobiomeStat")

  fx <- ap_da_planted(n_features = 30L, per = 12L)
  counts <- SummarizedExperiment::assay(fx$x, "counts")
  hashy <- paste0(sample(c(0:9, letters[1:6]), nrow(counts), replace = TRUE),
                  sprintf("%031d", seq_len(nrow(counts))))
  rownames(counts) <- hashy
  meta <- as.data.frame(SummarizedExperiment::colData(fx$x))
  y <- ap_import(counts, meta)

  da <- suppressMessages(ap_da(y, "dx", method = c("maaslin2", "linda"),
                               reference = "control", seed = 1L))
  for (m in unique(da$results$method)) {
    got <- da$results$feature[da$results$method == m]
    expect_setequal(got, hashy)
  }
  # Which is what makes the union equal the tested set rather than double it.
  expect_equal(length(unique(da$results$feature)), da$n_features)
})
