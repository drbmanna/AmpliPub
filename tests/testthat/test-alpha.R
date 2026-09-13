test_that("Hill numbers on a known community equal the answer worked out by hand", {
  # Four equally abundant species: every Hill number is 4, by definition,
  # because q only changes how abundance is weighted and here it is uniform.
  even <- matrix(c(25, 25, 25, 25), ncol = 1,
                 dimnames = list(paste0("f", 1:4), "s1"))
  res <- ap_alpha_one(even, c("q0", "q1", "q2", "evenness"), NULL)
  v <- stats::setNames(res$value, res$metric)
  expect_equal(unname(v["q0"]), 4)
  expect_equal(unname(v["q1"]), 4)
  expect_equal(unname(v["q2"]), 4)
  expect_equal(unname(v["evenness"]), 1)
})

test_that("Hill numbers order q0 >= q1 >= q2 on an uneven community", {
  uneven <- matrix(c(97, 1, 1, 1), ncol = 1,
                   dimnames = list(paste0("f", 1:4), "s1"))
  res <- ap_alpha_one(uneven, c("q0", "q1", "q2"), NULL)
  v <- stats::setNames(res$value, res$metric)
  expect_equal(unname(v["q0"]), 4)
  expect_gt(v["q0"], v["q1"])
  expect_gt(v["q1"], v["q2"])
  # Dominated by one species, so the effective number is near 1.
  expect_lt(v["q2"], 1.3)
})

test_that("Shannon entropy is reported in the base asked for", {
  m <- matrix(c(25, 25, 25, 25), ncol = 1, dimnames = list(paste0("f", 1:4), "s1"))
  expect_equal(unname(ap_shannon_entropy(m, base = 2)), 2)      # 4 species = 2 bits
  expect_equal(unname(ap_shannon_entropy(m, base = exp(1))), log(4))
})

test_that("a zero-abundance feature does not make entropy NaN", {
  m <- matrix(c(50, 50, 0), ncol = 1, dimnames = list(paste0("f", 1:3), "s1"))
  expect_equal(unname(ap_shannon_entropy(m, base = 2)), 1)
})

test_that("evenness is NA rather than Inf when a sample holds one feature", {
  m <- matrix(c(100, 0, 0), ncol = 1, dimnames = list(paste0("f", 1:3), "s1"))
  res <- ap_alpha_one(m, "evenness", NULL)
  expect_true(is.na(res$value))
})

# --- rarefaction ---

test_that("rarefaction draws exactly the requested depth from every sample", {
  counts <- ap_fixture_counts(n_features = 20L, n_per_group = 4L, depth = 3000)
  set.seed(3)
  sub <- ap_rarefy_matrix(counts, 500)
  expect_true(all(colSums(sub) == 500))
  expect_true(all(sub <= counts))
})

test_that("rarefaction lowers richness and is reproducible under a seed", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a1 <- ap_alpha(x, metrics = "q0", depth = 500, n_iter = 5L, seed = 99L)
  a2 <- ap_alpha(x, metrics = "q0", depth = 500, n_iter = 5L, seed = 99L)
  expect_equal(a1$values$value, a2$values$value)

  raw <- ap_alpha(x, metrics = "q0", rarefy = FALSE)
  expect_true(all(a1$values$value <= raw$values$value))
})

test_that("samples below the depth are dropped, loudly", {
  counts <- ap_fixture_counts(n_features = 20L, n_per_group = 4L, depth = 2000)
  counts[, 1] <- round(counts[, 1] / 10)
  meta <- ap_fixture_metadata(counts)
  x <- ap_import(counts, meta)
  expect_warning(a <- ap_alpha(x, metrics = "q0", depth = 1500, n_iter = 2L),
                 "below depth")
  expect_equal(a$dropped, "S01")
  expect_false("S01" %in% a$values$sample_id)
})

test_that("a depth above every library is refused rather than returning nothing", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_error(ap_alpha(x, depth = 1e9, n_iter = 2L), "above every library size")
})

test_that("rarefying a non-integer table is refused", {
  counts <- ap_fixture_counts(n_features = 10L, n_per_group = 3L)
  rel <- sweep(counts, 2, colSums(counts), "/")
  meta <- ap_fixture_metadata(counts)
  x <- suppressWarnings(ap_import(rel, meta, expect_counts = FALSE))
  expect_error(ap_alpha(x, metrics = "q0"), "needs integer counts")
})

test_that("not rarefying unequal depths warns that q0 tracks depth", {
  counts <- ap_fixture_counts(n_features = 20L, n_per_group = 4L)
  counts[, 1] <- round(counts[, 1] / 4)
  x <- ap_import(counts, ap_fixture_metadata(counts))
  expect_warning(ap_alpha(x, metrics = "q0", rarefy = FALSE), "rises with depth")
})

test_that("an unknown metric is refused by name", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_error(ap_alpha(x, metrics = "shannon"), "shannon")
})

test_that("Faith's PD is dropped with a warning when there is no tree", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_warning(a <- ap_alpha(x, metrics = c("q0", "faith_pd"), n_iter = 2L),
                 "needs a phylogeny")
  expect_false("faith_pd" %in% a$values$metric)
})

# --- Good's coverage ---

test_that("Good's coverage refuses on singleton-free data and explains why", {
  counts <- ap_fixture_counts(n_features = 15L, n_per_group = 4L, depth = 2000)
  counts[counts == 1] <- 2
  x <- ap_import(counts, ap_fixture_metadata(counts))
  expect_error(ap_goods_coverage(x), "DADA2 and Deblur remove singletons")
  expect_equal(unname(ap_goods_coverage(x, force = TRUE)), rep(1, ncol(counts)))
})

test_that("Good's coverage computes when singletons are present", {
  m <- matrix(c(1, 1, 98, 0, 0, 100), ncol = 2,
              dimnames = list(paste0("f", 1:3), c("s1", "s2")))
  x <- ap_import(m, data.frame(g = c("a", "b"), row.names = c("s1", "s2")))
  cov <- ap_goods_coverage(x)
  expect_equal(unname(cov["s1"]), 1 - 2 / 100)
  expect_equal(unname(cov["s2"]), 1)
})

# --- Chao1 and ACE ---

# Unequal depths and a long rare tail, so singletons are real rather than
# produced by subsampling.
ap_singleton_table <- function() {
  set.seed(8)
  m <- vapply(c(800, 1200, 1500, 2000), function(d) {
    stats::rmultinom(1, d, stats::rgamma(60, shape = 0.3, rate = 0.01) + 0.01)[, 1]
  }, numeric(60))
  dimnames(m) <- list(sprintf("f%02d", 1:60), sprintf("s%d", 1:4))
  m[rowSums(m) > 0, ]
}
ap_singleton_meta <- function() {
  data.frame(g = c("a", "a", "b", "b"), row.names = sprintf("s%d", 1:4))
}

test_that("Chao1 reaches ap_alpha from vegan on the right sample, checked by hand", {
  # The estimate itself is vegan's. What is ours is the wiring: transposing the
  # table and picking the right row and column. Sample s1 has counts
  # 1, 1, 1, 2, 5, so S_obs = 5, F1 = 3, F2 = 1 and Chao1 = 5 + 3 * 2 / (2 * 2) = 6.5.
  # Sample s2 holds no singletons, so Chao1 = S_obs = 3.
  m <- matrix(c(1, 1, 1, 2, 5, 0,
                0, 4, 0, 7, 0, 3), ncol = 2,
              dimnames = list(paste0("f", 1:6), c("s1", "s2")))
  x <- ap_import(m, data.frame(g = c("a", "b"), row.names = c("s1", "s2")))
  a <- suppressWarnings(ap_alpha(x, metrics = c("chao1", "ace"), rarefy = FALSE))
  v <- a$values[a$values$metric == "chao1", ]
  expect_equal(v$value[v$sample_id == "s1"], 6.5)
  expect_equal(v$value[v$sample_id == "s2"], 3)
  expect_true(all(a$values$metric[is.na(a$values$value)] == "ace"))
})

test_that("with no singletons Chao1 and ACE are refused, and when forced they equal q0", {
  counts <- ap_fixture_counts(n_features = 15L, n_per_group = 4L, depth = 2000)
  counts[counts == 1] <- 2
  x <- ap_import(counts, ap_fixture_metadata(counts))
  expect_error(ap_alpha(x, metrics = "chao1", rarefy = FALSE), "No feature has a count of 1")
  expect_error(ap_alpha(x, metrics = "ace", rarefy = FALSE), "DADA2 and Deblur")

  a <- suppressWarnings(ap_alpha(x, metrics = c("q0", "chao1", "ace"), rarefy = FALSE,
                                 force = TRUE))
  w <- split(a$values$value, a$values$metric)
  expect_equal(w$chao1, w$q0)
  expect_equal(w$ace, w$q0)
})

test_that("an equal-depth table is refused because rarefaction makes singletons", {
  set.seed(9)
  r <- ap_rarefy_matrix(ap_singleton_table(), 700)
  expect_gt(sum(r == 1), 0)
  x <- ap_import(r, ap_singleton_meta())
  expect_error(ap_alpha(x, metrics = "chao1", rarefy = FALSE), "rarefied table")
  expect_warning(a <- ap_alpha(x, metrics = "chao1", rarefy = FALSE, force = TRUE),
                 "equal-depth")
  expect_true(all(is.finite(a$values$value)))
})

test_that("real singletons at unequal depths give estimates at or above observed richness", {
  m <- ap_singleton_table()
  expect_gt(sum(m == 1), 0)
  x <- ap_import(m, ap_singleton_meta())
  a <- suppressWarnings(ap_alpha(x, metrics = c("q0", "chao1", "ace"), rarefy = FALSE))
  w <- split(a$values$value, a$values$metric)
  expect_true(all(w$chao1 >= w$q0))
  expect_true(all(w$ace >= w$q0, na.rm = TRUE))
  expect_gt(max(w$chao1 - w$q0), 0)
})

test_that("rarefying a singleton-free table is refused before the subsamples make singletons", {
  counts <- ap_fixture_counts(n_features = 30L, n_per_group = 4L, depth = 3000)
  counts[counts == 1] <- 2
  x <- ap_import(counts, ap_fixture_metadata(counts))
  # The subsamples would hold singletons the table never had...
  set.seed(1)
  expect_gt(sum(ap_rarefy_matrix(SummarizedExperiment::assay(x, "counts"), 500) == 1), 0)
  # ...and the guard looks before rarefying, so it still refuses.
  expect_error(ap_alpha(x, metrics = "chao1", depth = 500, n_iter = 2L),
               "No feature has a count of 1")
})

# --- depth candidates ---

test_that("depth candidates report the trade between samples kept and reads used", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  cand <- ap_depth_candidates(x, candidates = c(100, 5000, 1e6))
  expect_equal(cand$samples_retained[cand$depth == 100], ncol(x))
  expect_equal(cand$samples_retained[cand$depth == 1e6], 0L)
  expect_true(all(diff(cand$samples_retained) <= 0))
})

# --- Faith's PD implementation ---

test_that("PD on a hand-drawn tree equals the branch lengths summed by hand", {
  tr <- ape::read.tree(text = "((A:1,B:2):3,C:4);")
  ids <- c("A", "B", "C")
  idx <- ap_pd_index(tr, ids)

  # A alone: its own branch plus the internal branch to the root, 1 + 3 = 4.
  m <- matrix(c(1, 0, 0), ncol = 1, dimnames = list(ids, "s"))
  expect_equal(unname(ap_faith_pd(m, idx)), 4)

  # A and B: 1 + 2 + 3 = 6.
  m <- matrix(c(1, 1, 0), ncol = 1, dimnames = list(ids, "s"))
  expect_equal(unname(ap_faith_pd(m, idx)), 6)

  # Everything: the whole tree, 1 + 2 + 3 + 4 = 10.
  m <- matrix(c(1, 1, 1), ncol = 1, dimnames = list(ids, "s"))
  expect_equal(unname(ap_faith_pd(m, idx)), 10)
  expect_equal(unname(ap_faith_pd(m, idx)), idx$total)
})

test_that("PD counts a shared branch once, not once per descendant", {
  tr <- ape::read.tree(text = "((A:1,B:1):10,C:1);")
  idx <- ap_pd_index(tr, c("A", "B", "C"))
  m <- matrix(c(1, 1, 0), ncol = 1, dimnames = list(c("A", "B", "C"), "s"))
  # 1 + 1 + 10, not 1 + 10 + 1 + 10.
  expect_equal(unname(ap_faith_pd(m, idx)), 12)
})

test_that("PD ignores abundance and responds only to presence", {
  tr <- ape::read.tree(text = "((A:1,B:2):3,C:4);")
  idx <- ap_pd_index(tr, c("A", "B", "C"))
  a <- matrix(c(1, 5, 0), ncol = 1, dimnames = list(c("A", "B", "C"), "s"))
  b <- matrix(c(900, 7, 0), ncol = 1, dimnames = list(c("A", "B", "C"), "s"))
  expect_equal(ap_faith_pd(a, idx), ap_faith_pd(b, idx))
})

test_that("a tree without branch lengths is refused for PD", {
  tr <- ape::read.tree(text = "((A,B),C);")
  expect_error(ap_pd_index(tr, c("A", "B", "C")), "no branch lengths")
})

test_that("the fast PD index agrees with pruning the tree per sample", {
  # The slow, obviously correct implementation, on a tree big enough that an
  # indexing mistake would show.
  set.seed(5)
  tr <- ape::rtree(40)
  ids <- tr$tip.label
  idx <- ap_pd_index(tr, ids)

  counts <- matrix(stats::rbinom(40 * 6, 1, 0.4), nrow = 40,
                   dimnames = list(ids, paste0("s", 1:6)))
  counts[1, ] <- 1  # every sample keeps at least one tip

  slow <- vapply(seq_len(ncol(counts)), function(j) {
    present <- ids[counts[, j] > 0]
    if (length(present) == length(ids)) return(sum(tr$edge.length))
    pruned <- ape::keep.tip(tr, present)
    sum(pruned$edge.length)
  }, numeric(1))

  expect_equal(unname(ap_faith_pd(counts, idx)), slow, tolerance = 1e-10)
})
