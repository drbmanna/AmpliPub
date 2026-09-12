test_that("flat key: value documents parse, nested ones are refused", {
  parsed <- ap_parse_flat_yaml(c(
    "uuid: abc",
    "type: FeatureTable[Frequency]",
    "format: BIOMV210DirFmt",
    "data-size: 558.9 KiB"
  ))
  expect_equal(parsed$uuid, "abc")
  expect_equal(parsed$type, "FeatureTable[Frequency]")
  expect_equal(parsed[["data-size"]], "558.9 KiB")

  expect_error(
    ap_parse_flat_yaml(c("a: 1", "  b: 2")),
    "nesting"
  )
})

test_that("a distance matrix with a non-zero diagonal is refused", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("\tA\tB", "A\t0.5\t0.3", "B\t0.3\t0"), path)
  expect_error(ap_read_q2_distance(path), "non-zero diagonal")
})

test_that("a non-square distance matrix is refused", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("\tA\tB\tC", "A\t0\t0.3\t0.4", "B\t0.3\t0\t0.2"), path)
  expect_error(ap_read_q2_distance(path), "must be square")
})

test_that("a valid distance matrix becomes a dist", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("\tA\tB\tC", "A\t0\t0.3\t0.4", "B\t0.3\t0\t0.2", "C\t0.4\t0.2\t0"), path)
  d <- ap_read_q2_distance(path)
  expect_s3_class(d, "dist")
  expect_equal(attr(d, "Size"), 3L)
  expect_equal(as.numeric(d), c(0.3, 0.4, 0.2))
})

test_that("ordination files parse into coordinates, eigenvalues and proportions", {
  path <- withr::local_tempfile(fileext = ".txt")
  writeLines(c(
    "Eigvals\t2",
    "4.0\t1.0",
    "",
    "Proportion explained\t2",
    "0.8\t0.2",
    "",
    "Species\t0\t0",
    "",
    "Site\t3\t2",
    "A\t0.1\t0.2",
    "B\t-0.3\t0.4",
    "C\t0.2\t-0.6",
    "",
    "Biplot\t0\t0",
    "",
    "Site constraints\t0\t0"
  ), path)

  ord <- ap_read_q2_ordination(path)
  expect_equal(ord$eigvals, c(4, 1))
  expect_equal(ord$prop_explained, c(0.8, 0.2))
  expect_equal(dim(ord$vectors), c(3L, 2L))
  expect_equal(rownames(ord$vectors), c("A", "B", "C"))
  expect_equal(colnames(ord$vectors), c("PC1", "PC2"))
  expect_equal(unname(ord$vectors["B", ]), c(-0.3, 0.4))
})

test_that("an ordination with more axes than eigenvalues is refused", {
  path <- withr::local_tempfile(fileext = ".txt")
  writeLines(c(
    "Eigvals\t1",
    "4.0",
    "",
    "Site\t2\t2",
    "A\t0.1\t0.2",
    "B\t-0.3\t0.4"
  ), path)
  expect_error(ap_read_q2_ordination(path), "internally inconsistent")
})

test_that("duplicate feature IDs in a taxonomy file are refused", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("Feature ID\tTaxon", "A\td__Bacteria", "A\td__Archaea"), path)
  expect_error(ap_read_q2_taxonomy(path), "duplicate feature IDs")
})
