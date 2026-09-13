# --- UniFrac (weighted from mia/rbiom, unweighted from AmpliPub) ---

test_that("UniFrac matches scikit-bio on pairs that do not span the root", {
  # Reference values from scikit-bio 0.6.2, QIIME 2's engine, run in the
  # qiime2-amplicon-2025.7 environment on 2026-09-13. rbiom 2.2.1 gave a
  # different unweighted value on every one of these pairs, which is why
  # unweighted UniFrac is not taken from it. The weighted values agreed.
  t1 <- ape::read.tree(text = "((A:1,B:1):1,C:3);")
  t2 <- ape::read.tree(text = "(((A:1,B:1):1,C:1):1,D:5);")
  x1 <- c("A", "B", "C")
  x2 <- c("A", "B", "C", "D")
  pair <- function(a, b, ids) matrix(c(a, b), ncol = 2, dimnames = list(ids, c("s1", "s2")))
  u <- function(m, tr) as.numeric(ap_unifrac(m, tr, weighted = FALSE))
  w <- function(m, tr) as.numeric(ap_unifrac(m, tr, weighted = TRUE))

  expect_equal(u(pair(c(1, 0, 0), c(1, 1, 0), x1), t1), 1 / 3)
  expect_equal(u(pair(c(1, 0, 0), c(0, 1, 0), x1), t1), 2 / 3)
  expect_equal(u(pair(c(1, 0, 0, 0), c(0, 1, 0, 0), x2), t2), 1 / 2)
  expect_equal(u(pair(c(1, 0, 0, 0), c(1, 0, 1, 0), x2), t2), 1 / 4)

  expect_equal(w(pair(c(1, 0, 0), c(1, 1, 0), x1), t1), 1)
  expect_equal(w(pair(c(2, 1, 0), c(1, 2, 0), x1), t1), 2 / 3)
  expect_equal(w(pair(c(1, 0, 0, 0), c(0, 1, 0, 0), x2), t2), 2)
  expect_equal(w(pair(c(1, 0, 0, 0), c(1, 0, 1, 0), x2), t2), 1.5)
})

test_that("unweighted UniFrac on a hand-drawn tree equals the fraction computed by hand", {
  tr <- ape::read.tree(text = "((A:1,B:1):1,C:3);")
  ids <- c("A", "B", "C")

  # s1 has A only (branches A=1, internal=1, total 2).
  # s2 has C only (branch C=3, total 3).
  # Shared branch length: 0. Union: 5. Unique/union = 5/5 = 1.
  m <- matrix(c(1, 0, 0,
                0, 0, 1), ncol = 2, dimnames = list(ids, c("s1", "s2")))
  expect_equal(as.numeric(ap_unifrac(m, tr, weighted = FALSE)), 1)

  # s1 = A, s2 = A and B. Shared 2 (A + internal), union 3 (A + B + internal).
  # Unique = 1, so distance = 1/3.
  m <- matrix(c(1, 0, 0,
                1, 1, 0), ncol = 2, dimnames = list(ids, c("s1", "s2")))
  expect_equal(as.numeric(ap_unifrac(m, tr, weighted = FALSE)), 1 / 3)
})

test_that("a sample compared with itself is at distance zero", {
  tr <- ape::read.tree(text = "((A:1,B:2):3,C:4);")
  m <- matrix(c(5, 3, 0, 5, 3, 0), ncol = 2,
              dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  expect_equal(as.numeric(ap_unifrac(m, tr, weighted = FALSE)), 0)
  expect_equal(as.numeric(ap_unifrac(m, tr, weighted = TRUE)), 0)
})

test_that("unweighted UniFrac ignores abundance and weighted UniFrac does not", {
  tr <- ape::read.tree(text = "((A:1,B:1):1,C:3);")
  # Same taxa present in both, wildly different proportions.
  m <- matrix(c(99, 1, 0,
                1, 99, 0), ncol = 2, dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  expect_equal(as.numeric(ap_unifrac(m, tr, weighted = FALSE)), 0)
  expect_gt(as.numeric(ap_unifrac(m, tr, weighted = TRUE)), 0)
})

test_that("weighted UniFrac is the raw branch-weighted difference computed by hand", {
  tr <- ape::read.tree(text = "((A:1,B:1):1,C:3);")
  m <- matrix(c(1, 0, 0,
                0, 1, 0), ncol = 2, dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  # Edge A: |1 - 0| * 1 = 1. Edge B: |0 - 1| * 1 = 1. Internal: |1 - 1| * 1 = 0.
  # Edge C: 0. Total 2. A normalised form would divide this down to at most 1.
  expect_equal(as.numeric(ap_unifrac(m, tr, weighted = TRUE)), 2)
})

test_that("UniFrac on a tree rerooted with ape::root() matches its Newick-read twin", {
  # mia's Faith's PD crashed on this kind of tree. UniFrac does not crash, and
  # this checks it does not quietly give a different answer either.
  counts <- ap_fixture_counts()
  tr <- ap_fixture_tree(counts)
  twin <- ape::read.tree(text = ape::write.tree(tr, digits = 17))
  for (w in c(FALSE, TRUE)) {
    expect_equal(as.numeric(ap_unifrac(counts, tr, weighted = w)),
                 as.numeric(ap_unifrac(counts, twin, weighted = w)),
                 tolerance = 1e-10)
  }
})

test_that("sample order in the distance follows the table, not the package", {
  tr <- ape::read.tree(text = "((A:1,B:1):1,C:3);")
  m <- matrix(c(1, 0, 0,
                1, 1, 0,
                0, 0, 1), ncol = 3,
              dimnames = list(c("A", "B", "C"), c("zeta", "alpha", "mid")))
  d <- ap_unifrac(m, tr, weighted = FALSE)
  expect_equal(attr(d, "Labels"), c("zeta", "alpha", "mid"))
  expect_equal(as.matrix(d)["zeta", "alpha"], 1 / 3)
})

test_that("a tree without branch lengths is refused for UniFrac", {
  tr <- ape::read.tree(text = "((A,B),C);")
  m <- matrix(c(1, 0, 0, 0, 0, 1), ncol = 2, dimnames = list(c("A", "B", "C"), c("s1", "s2")))
  expect_error(ap_unifrac(m, tr, weighted = FALSE), "no branch lengths")
})

# --- ap_beta ---

test_that("ap_beta computes every requested metric", {
  x <- ap_fixture_object()
  b <- ap_beta(x, metrics = c("bray_curtis", "jaccard", "unweighted_unifrac",
                              "weighted_unifrac", "aitchison"),
               rarefy = FALSE)
  expect_length(b$distances, 5L)
  expect_true(all(vapply(b$distances, inherits, logical(1), "dist")))
  expect_true(all(vapply(b$distances, function(d) attr(d, "Size"), integer(1)) == ncol(x)))
})

test_that("UniFrac metrics are dropped with a warning when there is no tree", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_warning(b <- ap_beta(x, rarefy = FALSE), "need a phylogeny")
  expect_false(any(grepl("unifrac", names(b$distances))))
})

test_that("an unknown metric is refused by name", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  expect_error(ap_beta(x, metrics = "manhattan"), "manhattan")
})

test_that("rarefaction in ap_beta is reproducible under a seed", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b1 <- ap_beta(x, metrics = "bray_curtis", depth = 1000, seed = 7L)
  b2 <- ap_beta(x, metrics = "bray_curtis", depth = 1000, seed = 7L)
  expect_equal(as.numeric(b1$distances$bray_curtis),
               as.numeric(b2$distances$bray_curtis))
})

test_that("Aitchison distance is unchanged when a sample is scaled up", {
  # The point of a compositional metric: doubling a sample's read depth changes
  # nothing about its composition, so the distance must not move.
  m <- matrix(c(10, 20, 30, 5, 50, 5), ncol = 2,
              dimnames = list(paste0("f", 1:3), c("s1", "s2")))
  d1 <- ap_aitchison(m)
  m2 <- m
  m2[, 1] <- m2[, 1] * 10
  d2 <- ap_aitchison(m2)
  expect_equal(as.numeric(d1), as.numeric(d2))
})

test_that("a pseudocount of zero is refused, because the CLR is undefined there", {
  m <- matrix(c(1, 0, 2, 3), ncol = 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  expect_error(ap_aitchison(m, pseudocount = 0), "must be positive")
})

test_that("a table that is almost all zeros warns that the pseudocount dominates", {
  set.seed(2)
  m <- matrix(0, nrow = 50, ncol = 6,
              dimnames = list(paste0("f", 1:50), paste0("s", 1:6)))
  m[cbind(sample(50, 6), 1:6)] <- 100
  expect_warning(ap_aitchison(m), "dominated by the pseudocount")
})

# --- ordination ---

test_that("PCoA reports variance explained that sums to 1 over positive eigenvalues", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  ord <- ap_ordinate(b, "bray_curtis", method = "pcoa")
  expect_equal(sum(ord$prop_explained[ord$eig > 0]), 1)
  expect_true(all(diff(ord$prop_explained[ord$eig > 0]) <= 1e-12))
  expect_equal(nrow(ord$coords), ncol(x))
})

test_that("PCoA on a Euclidean distance has no negative eigenvalue mass", {
  set.seed(3)
  pts <- matrix(stats::rnorm(60), ncol = 3,
                dimnames = list(paste0("s", 1:20), NULL))
  ord <- ap_ordinate(stats::dist(pts), method = "pcoa")
  expect_lt(ord$negative_eigenvalue_fraction, 1e-8)
})

test_that("PCoA recovers a planted two-cluster structure on the first axis", {
  set.seed(9)
  pts <- rbind(
    matrix(stats::rnorm(30, -5), ncol = 3),
    matrix(stats::rnorm(30, 5), ncol = 3)
  )
  rownames(pts) <- paste0("s", 1:20)
  ord <- ap_ordinate(stats::dist(pts), method = "pcoa")
  grp <- rep(c("a", "b"), each = 10)
  expect_gt(abs(stats::cor(ord$coords[, 1], as.numeric(factor(grp)))), 0.95)
  expect_gt(ord$prop_explained[1], 0.8)
})

test_that("axis labels carry the variance explained for PCoA and not for NMDS", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  pcoa <- ap_ordinate(b, "bray_curtis", method = "pcoa")
  expect_match(ap_axis_label(pcoa, 1), "PCo1 \\([0-9.]+%\\)")

  nmds <- suppressWarnings(ap_ordinate(b, "bray_curtis", method = "nmds", trymax = 5L))
  expect_equal(ap_axis_label(nmds, 1), "NMDS1")
  expect_false(is.na(nmds$stress))
})

test_that("asking for a metric the object does not hold is refused", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  expect_error(ap_ordinate(b, "jaccard"), "not in this ap_beta object")
})
