test_that("variables are classified into roles", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  scan <- ap_scan_metadata(x)
  roles <- stats::setNames(scan$variables$role, scan$variables$variable)

  expect_equal(unname(roles["group"]), "binary")
  expect_equal(unname(roles["age"]), "numeric")
  expect_equal(unname(roles["constant"]), "constant")
  expect_equal(scan$n_samples, 24L)
})

test_that("a run or batch variable is flagged by name", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  scan <- ap_scan_metadata(x)
  expect_true("batch_run" %in% scan$batch_candidates)
  expect_false("age" %in% scan$batch_candidates)
})

test_that("Site is treated as a batch proxy, because in multi-centre studies it is the only one", {
  meta <- data.frame(
    Site = rep(c("U Michigan", "Toronto"), each = 6),
    dx = rep(c("normal", "cancer"), 6),
    row.names = sprintf("S%02d", 1:12),
    stringsAsFactors = FALSE
  )
  scan <- ap_scan_metadata(meta)
  expect_true("Site" %in% scan$batch_candidates)
})

test_that("repeated subject IDs are flagged as possible repeated measures", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  scan <- ap_scan_metadata(x, max_levels = 5L)
  expect_true("subject" %in% scan$repeated_measures)
  expect_false("group" %in% scan$repeated_measures)
})

test_that("a variable that tracks the group perfectly is screened as a confounder", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  scan <- ap_scan_metadata(x, group = "group")
  conf <- scan$confounders
  expect_true("batch_run" %in% conf$variable)
  expect_lt(conf$p_adj[conf$variable == "batch_run"], 0.05)
  # sex alternates independently of group and must not be flagged.
  expect_gt(conf$p_adj[conf$variable == "sex"], 0.05)
})

test_that("naming a group that does not exist is refused with the options listed", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_error(ap_scan_metadata(x, group = "nope"), "not in the metadata")
})

test_that("the scan prints without error", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  scan <- ap_scan_metadata(x, group = "group")
  expect_message(print(scan), "Metadata assessment")
})

# --- the grouping-variable guard ----------------------------------------------
#
# The failure this guard exists for is silent: name a column that cannot group
# samples and every stage still runs, the report still renders, and the
# interpretation sentences attach real numbers to a meaningless comparison.
# Each refusal below is therefore asserted to actually fire.

# cli wraps long messages, so a fragment can be split across lines. Normalise
# whitespace before matching rather than testing the wrapping.
ap_err_text <- function(expr) {
  e <- tryCatch(expr, error = function(e) e)
  if (!inherits(e, "condition")) return("")
  gsub("[[:space:]]+", " ", paste(conditionMessage(e), collapse = " "))
}

ap_guard_meta <- function() {
  data.frame(
    dx = rep(c("case", "control"), each = 6),
    sample_id = sprintf("S%02d", 1:12),
    study = "baxter",
    row.names = sprintf("S%02d", 1:12),
    stringsAsFactors = FALSE
  )
}

test_that("a constant grouping variable is refused, not silently analysed", {
  txt <- ap_err_text(ap_scan_metadata(ap_guard_meta(), group = "study"))
  expect_true(grepl("one value across every sample", txt, fixed = TRUE))
  expect_true(grepl("needs two groups", txt, fixed = TRUE))
})

test_that("an identifier-like grouping variable is refused", {
  txt <- ap_err_text(ap_scan_metadata(ap_guard_meta(), group = "sample_id"))
  expect_true(grepl("about one per sample", txt, fixed = TRUE))
  expect_true(grepl("sample or subject identifier", txt, fixed = TRUE))
})

test_that("the refusal names the variables that could be used instead", {
  txt <- ap_err_text(ap_scan_metadata(ap_guard_meta(), group = "study"))
  expect_true(grepl("Variables that can group these samples: dx", txt, fixed = TRUE))
})

test_that("force = TRUE downgrades the refusal to a warning and records it", {
  expect_warning(
    scan <- ap_scan_metadata(ap_guard_meta(), group = "study", force = TRUE),
    "force"
  )
  expect_equal(scan$group_status, "constant")
  expect_true(any(grepl("Forced", scan$group_notes)))
})

test_that("a usable grouping variable passes clean", {
  scan <- ap_scan_metadata(ap_guard_meta(), group = "dx")
  expect_equal(scan$group_status, "ok")
  expect_length(scan$group_notes, 0L)
})

test_that("a level smaller than three warns but does not stop the run", {
  meta <- data.frame(
    dx = rep(c("case", "control"), c(10, 2)),
    row.names = sprintf("S%02d", 1:12),
    stringsAsFactors = FALSE
  )
  expect_warning(scan <- ap_scan_metadata(meta, group = "dx"), "smallest level")
  expect_equal(scan$group_status, "warned")
  expect_true(any(grepl("dispersion test is degenerate", scan$group_notes)))
})

test_that("more levels than max_levels warns but does not stop the run", {
  meta <- data.frame(
    site = rep(sprintf("site%d", 1:6), each = 2),
    row.names = sprintf("S%02d", 1:12),
    stringsAsFactors = FALSE
  )
  expect_warning(
    scan <- ap_scan_metadata(meta, group = "site", max_levels = 5L),
    "6 levels"
  )
  expect_equal(scan$group_status, "warned")
})

test_that("guard notes reach the printed scan, which is what the report captures", {
  meta <- data.frame(
    dx = rep(c("case", "control"), c(10, 2)),
    row.names = sprintf("S%02d", 1:12),
    stringsAsFactors = FALSE
  )
  suppressWarnings(scan <- ap_scan_metadata(meta, group = "dx"))
  # cli prints to stderr, which is why the report captures the message stream.
  out <- paste(utils::capture.output(print(scan), type = "message"), collapse = " ")
  expect_true(grepl("Grouping variable", out, fixed = TRUE))
  expect_true(grepl("dispersion", out, fixed = TRUE))
})
