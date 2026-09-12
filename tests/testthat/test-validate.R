# Every guard gets a test that watches it fire. A guard nobody has seen fail is
# a guard we do not know works.

test_that("a non-matrix table is refused", {
  expect_error(ap_check_count_matrix(list(1, 2)), "must be a matrix")
})

test_that("a character matrix is refused", {
  m <- matrix("1", 2, 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_error(ap_check_count_matrix(m), "must be numeric")
})

test_that("missing values are refused, because an absent taxon is a zero", {
  m <- matrix(c(1, NA, 3, 4), 2, 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_error(ap_check_count_matrix(m), "missing value")
})

test_that("negative counts are refused", {
  m <- matrix(c(1, -2, 3, 4), 2, 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_error(ap_check_count_matrix(m), "negative")
})

test_that("a table with no IDs is refused", {
  m <- matrix(1:4, 2, 2)
  expect_error(ap_check_count_matrix(m), "no feature IDs")
})

test_that("duplicate feature IDs are refused", {
  m <- matrix(1:4, 2, 2, dimnames = list(c("a", "a"), c("s1", "s2")))
  expect_error(ap_check_count_matrix(m), "duplicate feature IDs")
})

test_that("duplicate sample IDs are refused", {
  m <- matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("s1", "s1")))
  expect_error(ap_check_count_matrix(m), "duplicate sample IDs")
})

test_that("a relative abundance table is recognised, not silently accepted", {
  m <- matrix(c(0.25, 0.75, 0.4, 0.6), 2, 2,
              dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_warning(ap_check_counts_are_integers(m), "relative abundance")
  expect_error(ap_check_counts_are_integers(m, strict = TRUE), "relative abundance")
})

test_that("integer counts pass quietly", {
  m <- matrix(c(3, 7, 4, 6), 2, 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_silent(ap_check_counts_are_integers(m))
})

test_that("a sample with zero total counts is refused", {
  m <- matrix(c(3, 7, 0, 0), 2, 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_error(ap_check_no_empty(m), "zero total counts")
})

test_that("features that are zero everywhere are dropped and the drop is announced", {
  m <- matrix(c(3, 0, 4, 0), 2, 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_message(out <- ap_check_no_empty(m), "Dropping 1 feature")
  expect_equal(rownames(out), "a")
})

# --- the join, the guard that matters most ---

test_that("a join with no shared IDs is refused and says what it saw", {
  expect_error(
    ap_check_join(c("s1", "s2"), c("x1", "x2")),
    "No sample IDs are shared"
  )
})

test_that("dropping more than max_drop of the table fails rather than proceeding", {
  table_ids <- sprintf("s%02d", 1:10)
  meta_ids <- sprintf("s%02d", 1:5)
  expect_error(
    ap_check_join(table_ids, meta_ids, max_drop = 0.1),
    "above the 10% limit"
  )
})

test_that("a small drop warns and reports both sides", {
  table_ids <- sprintf("s%02d", 1:10)
  meta_ids <- c(sprintf("s%02d", 2:10), "extra1")
  expect_warning(
    expect_message(
      res <- ap_check_join(table_ids, meta_ids, max_drop = 0.2),
      "1 metadata sample"
    ),
    "1 sample in the table"
  )
  expect_equal(res$only_table, "s01")
  expect_equal(res$only_meta, "extra1")
  expect_length(res$shared, 9L)
})

test_that("a feature absent from the tree is refused, because UniFrac is undefined for it", {
  tr <- ape::read.tree(text = "((A:0.1,B:0.2):0.3,C:0.4);")
  expect_error(ap_check_tree_covers(tr, c("A", "B", "D")), "absent from the tree")
})

test_that("a tree covering every feature passes", {
  tr <- ape::read.tree(text = "((A:0.1,B:0.2):0.3,C:0.4);")
  expect_true(ap_check_tree_covers(tr, c("A", "B")))
})

test_that("ap_validate passes a well-formed object", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_silent(ap_validate(x))
})

test_that("ap_validate catches a missing value written into the assay after import", {
  # SummarizedExperiment enforces the colData/assay correspondence itself, but it
  # is happy to hold NA counts. This is the reachable way to corrupt an object
  # after import, so it is the one worth guarding.
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  counts <- SummarizedExperiment::assay(x, "counts")
  counts[1, 1] <- NA_real_
  SummarizedExperiment::assay(x, "counts") <- counts
  expect_error(ap_validate(x), "missing value")
})

test_that("ap_validate catches a negative count written into the assay after import", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  counts <- SummarizedExperiment::assay(x, "counts")
  counts[2, 3] <- -1
  SummarizedExperiment::assay(x, "counts") <- counts
  expect_error(ap_validate(x), "negative")
})
