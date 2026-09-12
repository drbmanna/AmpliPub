# Synthetic fixtures.
#
# Small, deterministic, and built so the right answer is known before any
# function runs. Real data lives outside this repo and is never committed.

# A table with a planted effect: `effect_features` are enriched in group b by
# `effect_size` fold. Everything else is drawn from the same distribution in
# both groups, so a differential abundance method that calls anything outside
# `effect_features` is producing a false positive we can count.
ap_fixture_counts <- function(n_features = 30L, n_per_group = 12L,
                              effect_features = 1:5, effect_size = 8,
                              depth = 5000, seed = 42) {
  set.seed(seed)
  n_samples <- 2L * n_per_group
  effect_features <- effect_features[effect_features <= n_features]
  stopifnot(length(effect_features) > 0L)
  # The +2 floor keeps every feature detectable at this depth. Without it the
  # gamma tail produces features that are zero in every sample, which the import
  # guard then correctly drops, and the fixture would no longer have the
  # dimensions the tests assert.
  base <- stats::rgamma(n_features, shape = 0.7, rate = 0.02) + 2
  names(base) <- sprintf("ASV%02d", seq_len(n_features))

  group <- rep(c("a", "b"), each = n_per_group)
  counts <- matrix(0L, nrow = n_features, ncol = n_samples,
                   dimnames = list(names(base), sprintf("S%02d", seq_len(n_samples))))

  for (j in seq_len(n_samples)) {
    lambda <- base
    if (group[j] == "b") lambda[effect_features] <- lambda[effect_features] * effect_size
    probs <- lambda / sum(lambda)
    counts[, j] <- stats::rmultinom(1, size = depth, prob = probs)[, 1]
  }

  attr(counts, "group") <- group
  attr(counts, "effect_features") <- names(base)[effect_features]
  counts
}

ap_fixture_metadata <- function(counts) {
  group <- attr(counts, "group")
  n <- ncol(counts)
  set.seed(7)
  data.frame(
    group = group,
    # Tracks group perfectly: a planted confounder the scan must flag.
    batch_run = ifelse(group == "a", "run1", "run2"),
    # Independent of group.
    sex = rep(c("f", "m"), length.out = n),
    age = round(stats::runif(n, 30, 70)),
    subject = rep(sprintf("subj%02d", seq_len(n / 2)), times = 2),
    constant = "same",
    row.names = colnames(counts),
    stringsAsFactors = FALSE
  )
}

ap_fixture_tree <- function(counts, seed = 11) {
  set.seed(seed)
  tr <- ape::rtree(nrow(counts), tip.label = rownames(counts))
  ape::root(tr, outgroup = rownames(counts)[1], resolve.root = TRUE)
}

ap_fixture_taxonomy <- function(counts) {
  ids <- rownames(counts)
  n <- length(ids)
  phyla <- rep(c("Bacteroidota", "Bacillota_A_368345"), length.out = n)
  genera <- sprintf("Genus%02d", seq_len(n))
  taxon <- sprintf("d__Bacteria; p__%s; c__Cls; o__Ord; f__Fam; g__%s", phyla, genera)
  # A third of features stop at family, the way a real classifier truncates.
  trunc <- seq(3, n, by = 3)
  taxon[trunc] <- sprintf("d__Bacteria; p__%s; c__Cls; o__Ord; f__Fam", phyla[trunc])
  data.frame(feature_id = ids, taxon = taxon,
             confidence = seq(0.7, 0.99, length.out = n),
             stringsAsFactors = FALSE)
}

ap_fixture_object <- function(tree = TRUE, taxonomy = TRUE, ...) {
  counts <- ap_fixture_counts(...)
  meta <- ap_fixture_metadata(counts)
  ap_import(
    table = counts,
    metadata = meta,
    tree = if (tree) ap_fixture_tree(counts) else NULL,
    taxonomy = if (taxonomy) ap_fixture_taxonomy(counts) else NULL
  )
}

# Builds a real .qza on disk: a zip whose single top-level directory is the
# UUID, holding metadata.yaml, VERSION and data/. Layout verified against
# QIIME 2 amplicon 2025.7.0, archive format 7.0.
ap_fixture_qza <- function(dir, type, format, files,
                           uuid = "11111111-2222-3333-4444-555555555555",
                           framework = "2025.7.0", archive = "7.0") {
  skip_if_not_installed("zip")
  root <- file.path(dir, uuid)
  dir.create(file.path(root, "data"), recursive = TRUE, showWarnings = FALSE)

  writeLines(
    c(paste0("uuid: ", uuid), paste0("type: ", type), paste0("format: ", format)),
    file.path(root, "metadata.yaml")
  )
  writeLines(
    c("QIIME 2", paste0("archive: ", archive), paste0("framework: ", framework)),
    file.path(root, "VERSION")
  )
  for (nm in names(files)) {
    writeLines(files[[nm]], file.path(root, "data", nm))
  }

  out <- file.path(dir, paste0(uuid, ".qza"))
  zip::zip(out, files = uuid, root = dir, mode = "cherry-pick")
  out
}
