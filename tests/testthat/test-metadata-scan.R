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
