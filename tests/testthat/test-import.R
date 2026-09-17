test_that("import from an in-memory matrix produces a valid object", {
  x <- ap_fixture_object()
  expect_s4_class(x, "TreeSummarizedExperiment")
  expect_equal(dim(x), c(30L, 24L))
  expect_true("counts" %in% SummarizedExperiment::assayNames(x))
  expect_silent(ap_validate(x))
})

test_that("the provenance log records what was read and what was done", {
  x <- ap_fixture_object()
  log <- ap_provenance(x)
  expect_setequal(vapply(log$inputs, `[[`, character(1), "role"),
                  c("table", "metadata", "taxonomy", "tree"))
  expect_equal(log$steps[[1]]$step, "import")
  expect_equal(log$steps[[1]]$detail$samples, 24L)
  expect_true(!is.null(log$session$r_version))
})

test_that("an object with no AmpliPub log says so rather than returning nothing", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1, 2, 2, dimnames = list(c("a", "b"), c("s1", "s2"))))
  )
  expect_error(ap_provenance(se), "no AmpliPub provenance log")
})

test_that("taxonomy is parsed into rank columns at import", {
  x <- ap_fixture_object()
  rd <- SummarizedExperiment::rowData(x)
  expect_true(all(c("taxon", "phylum", "genus", "confidence") %in% colnames(rd)))
  # Every third feature is truncated at family in the fixture.
  expect_equal(sum(is.na(rd$genus)), length(seq(3, nrow(x), by = 3)))
})

test_that("the tree is pruned to the features the table actually holds", {
  counts <- ap_fixture_counts()
  tr <- ap_fixture_tree(counts)
  smaller <- counts[1:10, ]
  meta <- ap_fixture_metadata(counts)
  expect_message(
    x <- ap_import(smaller, meta, tree = tr),
    "Pruning 20 tree tips"
  )
  expect_equal(length(TreeSummarizedExperiment::rowTree(x)$tip.label), 10L)
})

# --- sample IDs must survive the round trip as text ---

test_that("numeric-looking sample IDs are not reformatted by the metadata reader", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("sample-id\tgroup", "007\ta", "2003650\tb", "0012\ta"), path)
  meta <- ap_load_metadata(path)
  expect_equal(rownames(meta), c("007", "2003650", "0012"))
  expect_type(meta$group, "character")
})

test_that("a QIIME 2 #q2:types directive row is not read as a sample", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("sample-id\tgroup\tage",
               "#q2:types\tcategorical\tnumeric",
               "S1\ta\t30",
               "S2\tb\t40"), path)
  meta <- ap_load_metadata(path)
  expect_equal(rownames(meta), c("S1", "S2"))
  expect_type(meta$age, "double")
  expect_type(meta$group, "character")
})

test_that("column types are inferred per column, not guessed from the first value", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("sample-id\tmixed\tnum", "S1\t1\t1.5", "S2\tabc\t2.5"), path)
  meta <- ap_load_metadata(path)
  expect_type(meta$mixed, "character")
  expect_type(meta$num, "double")
})

test_that("duplicate sample IDs in a metadata file are refused", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("sample-id\tgroup", "S1\ta", "S1\tb"), path)
  expect_error(ap_load_metadata(path), "duplicate sample IDs")
})

test_that("empty strings become NA rather than a level called empty", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("sample-id\tgroup", "S1\ta", "S2\t"), path)
  meta <- ap_load_metadata(path)
  expect_true(is.na(meta$group[2]))
})

# --- table readers ---

test_that("a TSV feature table reads with the same values as the matrix it came from", {
  counts <- ap_fixture_counts(n_features = 5L, n_per_group = 3L)
  path <- withr::local_tempfile(fileext = ".tsv")
  df <- data.frame(`#OTU ID` = rownames(counts), counts, check.names = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)

  m <- ap_load_table(path)
  expect_equal(unname(m), unname(counts[, , drop = TRUE]), ignore_attr = TRUE)
  expect_equal(rownames(m), rownames(counts))
})

test_that("a biom-convert banner line above the header is skipped", {
  counts <- ap_fixture_counts(n_features = 4L, n_per_group = 2L)
  path <- withr::local_tempfile(fileext = ".tsv")
  con <- file(path, "wt")
  writeLines("# Constructed from biom file", con)
  df <- data.frame(`#OTU ID` = rownames(counts), counts, check.names = FALSE)
  utils::write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)

  m <- ap_load_table(path)
  expect_equal(dim(m), dim(counts))
  expect_equal(rownames(m), rownames(counts))
})

test_that("a taxonomy column left in the feature table is refused with the fix named", {
  path <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("#OTU ID\tS1\tS2\ttaxonomy",
               "ASV1\t5\t3\td__Bacteria",
               "ASV2\t2\t9\td__Archaea"), path)
  expect_error(ap_load_table(path), "taxonomy")
})

test_that("features_are_rows = FALSE transposes rather than failing", {
  counts <- ap_fixture_counts(n_features = 6L, n_per_group = 3L)
  m <- ap_load_table(t(counts), features_are_rows = FALSE)
  expect_equal(dim(m), dim(counts))
  expect_equal(rownames(m), rownames(counts))
})

test_that("import fails when too many samples lack metadata", {
  counts <- ap_fixture_counts(n_features = 10L, n_per_group = 5L)
  meta <- ap_fixture_metadata(counts)[1:3, , drop = FALSE]
  expect_error(ap_import(counts, meta), "above the 10% limit")
})


# ---- feature order is imposed, not inherited ------------------------------
#
# Upstream order is not stable: the same Snakemake config run twice produced the same 829
# features in a different order, and ALDEx2 draws its Monte Carlo instances feature by
# feature from the RNG stream, so an identical seed on a differently ordered table gave
# different results. A seed cannot fix that, so `ap_import()` imposes the order.

test_that("features come back in a fixed order whatever order they arrived in", {
  counts <- ap_fixture_counts()
  meta <- ap_fixture_metadata(counts)

  set.seed(99)
  shuffled <- counts[sample(nrow(counts)), , drop = FALSE]
  expect_false(identical(rownames(counts), rownames(shuffled)))

  a <- ap_import(counts, meta)
  b <- ap_import(shuffled, meta)

  expect_identical(rownames(a), rownames(b))
  expect_identical(rownames(a), sort(rownames(counts), method = "radix"))
  # the values must travel with their feature, not just the names
  expect_equal(SummarizedExperiment::assay(a, "counts"),
               SummarizedExperiment::assay(b, "counts"))
})

test_that("a shuffled table gives the same differential abundance result", {
  counts <- ap_fixture_counts()
  meta <- ap_fixture_metadata(counts)
  set.seed(7)
  shuffled <- counts[sample(nrow(counts)), , drop = FALSE]

  a <- suppressWarnings(ap_da(ap_import(counts, meta), group = "group",
                              method = "aldex2", seed = 1L))
  b <- suppressWarnings(ap_da(ap_import(shuffled, meta), group = "group",
                              method = "aldex2", seed = 1L))
  key <- function(d) d$results[order(d$results$feature), c("feature", "effect", "p")]
  expect_equal(key(a), key(b), tolerance = 1e-12)
})

test_that("taxonomy and tree follow the reordered features", {
  counts <- ap_fixture_counts()
  meta <- ap_fixture_metadata(counts)
  tax <- ap_fixture_taxonomy(counts)
  set.seed(11)
  shuffled <- counts[sample(nrow(counts)), , drop = FALSE]

  x <- ap_import(shuffled, meta, taxonomy = tax)
  rd <- as.data.frame(SummarizedExperiment::rowData(x))
  direct <- ap_import(counts, meta, taxonomy = tax)
  expect_identical(rownames(rd), rownames(as.data.frame(SummarizedExperiment::rowData(direct))))
  expect_identical(rd$genus,
                   as.data.frame(SummarizedExperiment::rowData(direct))$genus)
})
