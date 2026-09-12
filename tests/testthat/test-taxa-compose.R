test_that("collapsing preserves every read", {
  x <- ap_fixture_object()
  cx <- ap_collapse(x, "genus")
  expect_equal(sum(SummarizedExperiment::assay(cx, "counts")),
               sum(SummarizedExperiment::assay(x, "counts")))
  expect_equal(ncol(cx), ncol(x))
  expect_lte(nrow(cx), nrow(x))
})

test_that("collapsing to phylum yields exactly the phyla in the fixture", {
  x <- ap_fixture_object()
  cx <- ap_collapse(x, "phylum")
  expect_setequal(rownames(cx), c("p__Bacteroidota", "p__Bacillota_A_368345"))
  # Two phyla, so each sample's two values must sum to its original depth.
  expect_equal(colSums(SummarizedExperiment::assay(cx, "counts")),
               colSums(SummarizedExperiment::assay(x, "counts")))
})

test_that("features unassigned at the rank are kept under a label that says so", {
  x <- ap_fixture_object()
  cx <- ap_collapse(x, "genus")
  expect_true(any(grepl("unassigned genus", rownames(cx), fixed = TRUE)))
  expect_equal(sum(SummarizedExperiment::assay(cx, "counts")),
               sum(SummarizedExperiment::assay(x, "counts")))
})

test_that("dropping the unassigned is announced with the reads it costs", {
  x <- ap_fixture_object()
  expect_message(cx <- ap_collapse(x, "genus", drop_unassigned = TRUE),
                 "Every remaining proportion")
  expect_false(any(grepl("Unassigned", rownames(cx), fixed = TRUE)))
})

test_that("collapsing records the step in the provenance log", {
  x <- ap_fixture_object()
  cx <- ap_collapse(x, "genus")
  steps <- vapply(ap_provenance(cx)$steps, `[[`, character(1), "step")
  expect_true("collapse" %in% steps)
  detail <- ap_provenance(cx)$steps[[which(steps == "collapse")]]$detail
  expect_equal(detail$rank, "genus")
  expect_true(detail$reads_unchanged)
})

test_that("collapsing without taxonomy is refused, naming what is available", {
  x <- ap_fixture_object(taxonomy = FALSE)
  expect_error(ap_collapse(x, "genus"), "imported without")
})

# --- top taxa ---

test_that("every sample's proportions sum to 1 once Other is included", {
  x <- ap_fixture_object()
  df <- ap_top_taxa(ap_collapse(x, "genus"), n = 5L)
  sums <- tapply(df$relative_abundance, df$sample_id, sum)
  expect_true(all(abs(sums - 1) < 1e-10))
})

test_that("the pooled remainder is reported, not hidden", {
  x <- ap_fixture_object()
  cx <- ap_collapse(x, "genus")
  df <- ap_top_taxa(cx, n = 3L)
  expect_equal(attr(df, "n_kept"), 3L)
  expect_equal(attr(df, "n_pooled"), nrow(cx) - 3L)
  expect_gt(attr(df, "pooled_mean_abundance"), 0)
  expect_true("Other" %in% levels(df$taxon))
})

test_that("Other is the last level so it stacks apart from real taxa", {
  x <- ap_fixture_object()
  df <- ap_top_taxa(ap_collapse(x, "genus"), n = 4L)
  expect_equal(utils::tail(levels(df$taxon), 1), "Other")
})

test_that("taxa are ranked by mean relative abundance, not by total counts", {
  # One deeply sequenced sample carries a taxon that is rare everywhere else.
  # Ranking on totals would promote it; ranking on the mean proportion must not.
  m <- matrix(0, nrow = 3, ncol = 4,
              dimnames = list(c("common", "deep_only", "rare"), paste0("s", 1:4)))
  m["common", ] <- c(90, 90, 90, 900)
  m["deep_only", ] <- c(5, 5, 5, 9000)
  m["rare", ] <- c(5, 5, 5, 100)
  meta <- data.frame(g = rep("a", 4), row.names = colnames(m))
  x <- ap_import(m, meta)

  df <- ap_top_taxa(x, n = 1L)
  expect_equal(as.character(levels(df$taxon)[1]), "common")
})

test_that("a group argument keeps taxa that are abundant only in a small group", {
  m <- matrix(1, nrow = 4, ncol = 20,
              dimnames = list(paste0("t", 1:4), paste0("s", 1:20)))
  m["t1", ] <- 1000                     # abundant everywhere
  m["t2", ] <- 500                      # second everywhere
  m["t3", 1:3] <- 2000                  # abundant only in the small group
  meta <- data.frame(g = c(rep("small", 3), rep("big", 17)),
                     row.names = colnames(m), stringsAsFactors = FALSE)
  x <- ap_import(m, meta)

  without <- ap_top_taxa(x, n = 2L)
  with_group <- ap_top_taxa(x, n = 2L, group = "g")
  expect_false("t3" %in% levels(without$taxon))
  expect_true("t3" %in% levels(with_group$taxon))
})

test_that("asking for more taxa than exist keeps all of them and pools nothing", {
  x <- ap_fixture_object()
  cx <- ap_collapse(x, "phylum")
  df <- ap_top_taxa(cx, n = 100L)
  expect_equal(attr(df, "n_pooled"), 0L)
  expect_false("Other" %in% levels(df$taxon))
})

# --- plots ---

test_that("the stacked bar chart builds with one bar per sample", {
  x <- ap_fixture_object()
  p <- ap_plot_taxa_bar(x, rank = "genus", n = 5L, group = "group")
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("the group-mean variant builds and needs a group", {
  x <- ap_fixture_object()
  expect_s3_class(ap_plot_taxa_bar(x, rank = "phylum", n = 5L, group = "group",
                                   mode = "group"), "ggplot")
  expect_error(ap_plot_taxa_bar(x, rank = "phylum", n = 5L, mode = "group"),
               "needs a `group` variable")
})

test_that("the caption states the pooling even after the group averaging reshape", {
  x <- ap_fixture_object()
  p <- ap_plot_taxa_bar(x, rank = "genus", n = 3L, group = "group", mode = "group")
  # The caption is wrapped for the device, so newlines are removed before
  # matching rather than the assertions being weakened to fit around them.
  flat <- gsub("\n", " ", p$labels$caption, fixed = TRUE)
  expect_match(flat, "pooled as Other")
  expect_match(flat, "mean relative abundance")
})

test_that("Other is drawn in grey rather than taking a categorical colour", {
  # Checked on the rendered fills, not the scale object's internals, so this
  # keeps working across ggplot2 versions.
  m <- matrix(c(80, 15, 5, 4, 1, 70, 20, 5, 4, 1), ncol = 2,
              dimnames = list(paste0("t", 1:5), c("s1", "s2")))
  x <- ap_import(m, data.frame(g = c("a", "b"), row.names = c("s1", "s2")))
  p <- ap_plot_taxa_bar(ap_top_taxa(x, n = 2L), rank = "feature")
  fills <- unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
  expect_true("#D9D9D9" %in% fills)
  expect_equal(sum(fills == "#D9D9D9"), 1L)
})

test_that("the heatmap builds and excludes the pooled remainder", {
  x <- ap_fixture_object()
  p <- ap_plot_taxa_heatmap(x, rank = "genus", n = 4L, group = "group")
  expect_s3_class(p, "ggplot")
  expect_false("Other" %in% levels(p$data$taxon))
  expect_silent(ggplot2::ggplot_build(p))
})
