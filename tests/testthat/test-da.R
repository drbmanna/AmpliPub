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

# ANCOM-BC2's robust call (2026-09-18). ANCOMBC defines diff_robust as q < alpha AND
# passed_ss, and recommends it for the final call. On Baxter 145 of 146 calls failed the
# sensitivity analysis and were being counted. These tests mock the package output so the
# guard is proven to fire without depending on a real run producing a failure.
fake_ancombc2_out <- function(taxa, lvl = "b") {
  n <- length(taxa)
  out <- data.frame(taxon = taxa)
  out[[paste0("lfc_group", lvl)]] <- c(2, 1.5, rep(0.01, n - 2))
  out[[paste0("se_group", lvl)]] <- 0.2
  out[[paste0("W_group", lvl)]] <- out[[paste0("lfc_group", lvl)]] / 0.2
  out[[paste0("p_group", lvl)]] <- c(1e-6, 1e-6, rep(0.9, n - 2))
  out[[paste0("q_group", lvl)]] <- c(1e-5, 1e-5, rep(0.95, n - 2))
  out[[paste0("passed_ss_group", lvl)]] <- c(TRUE, FALSE, rep(TRUE, n - 2))
  out
}

test_that("an ANCOM-BC2 call that failed its sensitivity analysis is not significant", {
  skip_if_not_installed("ANCOMBC")
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  taxa <- rownames(x)
  local_mocked_bindings(ancombc2 = function(...) list(res = fake_ancombc2_out(taxa)),
                        .package = "ANCOMBC")
  da <- suppressMessages(ap_da(x, "group", method = "ancombc2", prv_cut = 0))
  r <- da$results
  robust <- r[r$feature == taxa[1], ]
  fragile <- r[r$feature == taxa[2], ]
  expect_true(robust$significant)
  expect_false(robust$failed_sensitivity)
  # Same q-value, failed the sensitivity analysis: kept, flagged, not counted.
  expect_lt(fragile$p_adj, da$alpha)
  expect_false(fragile$significant)
  expect_true(fragile$failed_sensitivity)
  expect_equal(fragile$note, "failed pseudocount sensitivity")
  expect_equal(sum(r$significant), 1L)
})

test_that("the concordance counts only ANCOM-BC2's robust calls", {
  skip_if_not_installed("ANCOMBC")
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  taxa <- rownames(x)
  local_mocked_bindings(ancombc2 = function(...) list(res = fake_ancombc2_out(taxa)),
                        .package = "ANCOMBC")
  da <- suppressMessages(ap_da(x, "group", method = "ancombc2", prv_cut = 0))
  cc <- ap_da_concordance(da)
  expect_equal(as.integer(cc$n_by_method[["ancombc2"]]), 1L)
})

test_that("ancombc2 output without a passed_ss column is refused, not trusted", {
  skip_if_not_installed("ANCOMBC")
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  taxa <- rownames(x)
  bad <- fake_ancombc2_out(taxa)
  bad$passed_ss_groupb <- NULL
  local_mocked_bindings(ancombc2 = function(...) list(res = bad), .package = "ANCOMBC")
  expect_error(suppressMessages(ap_da(x, "group", method = "ancombc2", prv_cut = 0)),
               "Every method failed")
})

test_that("methods without their own robustness check are unaffected", {
  r <- ap_da_row(feature = c("f1", "f2"), method = "linda", effect = c(1, -1),
                 effect_scale = "log2 fold change", se = 0.1, statistic = 10,
                 p = 1e-6, p_adj = 1e-5, contrast = "b")
  expect_true(all(is.na(r$passed_sensitivity)))
})

# ANCOM-BC2 structural zeros (2026-09-19). ancombc2 drops a taxon absent from a whole group
# from `res` and lists it only in `zero_ind`; AmpliPub adds it back as a directional call.
test_that("a structural zero comes back as a call with a direction and no p-value", {
  skip_if_not_installed("ANCOMBC")
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  taxa <- rownames(x)
  est <- fake_ancombc2_out(taxa[-(1:2)])
  zi <- data.frame(taxon = taxa, check.names = FALSE)
  zi[["structural_zero (group = a)"]] <- c(TRUE, FALSE, rep(FALSE, length(taxa) - 2))
  zi[["structural_zero (group = b)"]] <- c(FALSE, TRUE, rep(FALSE, length(taxa) - 2))
  local_mocked_bindings(ancombc2 = function(...) list(res = est, zero_ind = zi),
                        .package = "ANCOMBC")
  da <- suppressMessages(ap_da(x, "group", method = "ancombc2", prv_cut = 0))
  r <- da$results
  absent_a <- r[r$feature == taxa[1], ]
  absent_b <- r[r$feature == taxa[2], ]
  # Reference is a. Absent from a = higher in b (+1); absent from b = lower (-1).
  expect_true(absent_a$structural_zero && absent_a$significant)
  expect_equal(absent_a$direction, 1)
  expect_equal(absent_b$direction, -1)
  expect_true(is.na(absent_a$p_adj) && is.na(absent_a$effect))
  expect_match(absent_a$note, "absent from a")
  # No feature vanished: every input taxon has an ANCOM-BC2 row.
  expect_setequal(r$feature, taxa)
  cc <- ap_da_concordance(da)
  expect_equal(cc$features$n_methods[cc$features$feature == taxa[1]], 1L)
})

test_that("real ancombc2 declares a planted structural zero and AmpliPub keeps it", {
  skip_on_cran()
  skip_if_no_ancombc2()
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  cnt <- SummarizedExperiment::assay(x, "counts")
  g <- SummarizedExperiment::colData(x)$group
  target <- rownames(cnt)[which.max(rowSums(cnt))]
  cnt[target, g == "a"] <- 0
  SummarizedExperiment::assay(x, "counts") <- cnt
  da <- suppressWarnings(suppressMessages(
    ap_da(x, "group", method = "ancombc2", prv_cut = 0, reference = "a")))
  row <- da$results[da$results$feature == target, ]
  expect_equal(nrow(row), 1L)
  expect_true(row$structural_zero)
  expect_equal(row$direction, 1)
  expect_true(row$significant)
})
