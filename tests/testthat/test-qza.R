test_that("ap_qza_info reads the header of a real archive layout", {
  dir <- withr::local_tempdir()
  path <- ap_fixture_qza(
    dir,
    type = "FeatureData[Taxonomy]",
    format = "TSVTaxonomyDirectoryFormat",
    files = list("taxonomy.tsv" = c(
      "Feature ID\tTaxon\tConfidence",
      "ASV01\td__Bacteria; p__Bacteroidota\t0.99"
    ))
  )

  info <- ap_qza_info(path)
  expect_equal(info$type, "FeatureData[Taxonomy]")
  expect_equal(info$format, "TSVTaxonomyDirectoryFormat")
  expect_equal(info$framework_version, "2025.7.0")
  expect_equal(info$archive_version, "7.0")
  expect_equal(info$data_files, "taxonomy.tsv")
})

test_that("a taxonomy artifact round-trips to a data frame", {
  dir <- withr::local_tempdir()
  path <- ap_fixture_qza(
    dir,
    type = "FeatureData[Taxonomy]",
    format = "TSVTaxonomyDirectoryFormat",
    files = list("taxonomy.tsv" = c(
      "Feature ID\tTaxon\tConfidence",
      "ASV01\td__Bacteria; p__Bacteroidota\t0.99",
      "ASV02\td__Bacteria; p__Bacillota_A_368345\t0.80"
    ))
  )

  tax <- ap_read_qza(path)
  expect_equal(names(tax), c("feature_id", "taxon", "confidence"))
  expect_equal(nrow(tax), 2L)
  expect_type(tax$confidence, "double")
  expect_equal(attr(tax, "qza_type"), "FeatureData[Taxonomy]")
  expect_equal(attr(tax, "qza_framework"), "2025.7.0")
})

test_that("a phylogeny artifact round-trips to a phylo", {
  dir <- withr::local_tempdir()
  path <- ap_fixture_qza(
    dir,
    type = "Phylogeny[Rooted]",
    format = "NewickDirectoryFormat",
    files = list("tree.nwk" = "((A:0.1,B:0.2):0.3,C:0.4);")
  )
  tr <- ap_read_qza(path)
  expect_s3_class(tr, "phylo")
  expect_setequal(tr$tip.label, c("A", "B", "C"))
})

# --- guards ---

test_that("a file that is not a zip is refused as not an artifact", {
  path <- withr::local_tempfile(fileext = ".qza")
  writeLines("this is not a zip archive", path)
  expect_error(ap_qza_info(path), "not a QIIME 2 artifact|readable zip")
})

test_that("a missing file is refused", {
  expect_error(ap_qza_info("no/such/file.qza"), "not found")
})

test_that("an empty file is refused", {
  path <- withr::local_tempfile(fileext = ".qza")
  file.create(path)
  expect_error(ap_qza_info(path), "empty")
})

test_that("a zip with more than one top-level directory is refused", {
  skip_if_not_installed("zip")
  dir <- withr::local_tempdir()
  for (d in c("aaa", "bbb")) {
    dir.create(file.path(dir, d))
    writeLines("x", file.path(dir, d, "f.txt"))
  }
  path <- file.path(dir, "two.qza")
  zip::zip(path, files = c("aaa", "bbb"), root = dir, mode = "cherry-pick")
  expect_error(ap_qza_info(path), "exactly one top-level directory")
})

test_that("a UUID that disagrees with metadata.yaml is refused", {
  skip_if_not_installed("zip")
  dir <- withr::local_tempdir()
  uuid <- "aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
  root <- file.path(dir, uuid)
  dir.create(file.path(root, "data"), recursive = TRUE)
  writeLines(c("uuid: 99999999-0000-0000-0000-000000000000",
               "type: FeatureData[Taxonomy]",
               "format: TSVTaxonomyDirectoryFormat"),
             file.path(root, "metadata.yaml"))
  writeLines("Feature ID\tTaxon", file.path(root, "data", "taxonomy.tsv"))
  path <- file.path(dir, "bad.qza")
  zip::zip(path, files = uuid, root = dir, mode = "cherry-pick")

  expect_error(ap_qza_info(path), "UUID mismatch")
})

test_that("an unsupported semantic type names itself in the error", {
  dir <- withr::local_tempdir()
  path <- ap_fixture_qza(
    dir,
    type = "SampleData[SequencesWithQuality]",
    format = "SingleLanePerSampleSingleEndFastqDirFmt",
    files = list("MANIFEST" = "sample-id,filename,direction")
  )
  expect_error(ap_read_qza(path), "SampleData\\[SequencesWithQuality\\]")
})

test_that("an artifact of the wrong type is refused by ap_import", {
  dir <- withr::local_tempdir()
  path <- ap_fixture_qza(
    dir,
    type = "Phylogeny[Rooted]",
    format = "NewickDirectoryFormat",
    files = list("tree.nwk" = "(A:0.1,B:0.2);")
  )
  expect_error(
    ap_import(path, data.frame(g = "x", row.names = "A")),
    "not a FeatureTable"
  )
})
