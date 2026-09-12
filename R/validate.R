# Guards.
#
# These run at import and at every boundary that changes the sample or feature
# set. The failure mode they exist for is not a crash. It is a join that
# silently keeps 40 samples out of 490 and returns a confident PERMANOVA on the
# survivors. Each guard has a test that proves it fires; a guard nobody has
# watched fail is a guard we do not know works.

#' @keywords internal
ap_check_count_matrix <- function(m, name = "feature table") {
  ap_assert(is.matrix(m), "The {name} must be a matrix, not {class(m)[1]}.")
  ap_assert(is.numeric(m), "The {name} must be numeric, not {typeof(m)}.")
  ap_assert(nrow(m) > 0L && ncol(m) > 0L,
            "The {name} is empty: {nrow(m)} features x {ncol(m)} samples.")

  n_na <- sum(is.na(m))
  ap_assert(n_na == 0L,
            paste0("The {name} holds {n_na} missing value{?s}. Counts cannot be missing; ",
                   "an absent taxon is a zero. Fix the table before importing."))

  n_neg <- sum(m < 0)
  ap_assert(n_neg == 0L, "The {name} holds {n_neg} negative value{?s}.")

  ap_assert(!is.null(rownames(m)), "The {name} has no feature IDs (row names).")
  ap_assert(!is.null(colnames(m)), "The {name} has no sample IDs (column names).")

  dup_f <- rownames(m)[duplicated(rownames(m))]
  ap_assert(length(dup_f) == 0L,
            "The {name} has duplicate feature IDs: {paste(utils::head(dup_f, 5), collapse = ', ')}.")
  dup_s <- colnames(m)[duplicated(colnames(m))]
  ap_assert(length(dup_s) == 0L,
            "The {name} has duplicate sample IDs: {paste(utils::head(dup_s, 5), collapse = ', ')}.")

  invisible(TRUE)
}

# Counts vs. something already normalised. Methods downstream (ANCOM-BC2,
# ALDEx2, rarefaction) are defined on counts, and handing them a relative
# abundance table produces numbers rather than an error.
#' @keywords internal
ap_check_counts_are_integers <- function(m, strict = FALSE) {
  frac <- m[m != round(m)]
  if (length(frac) == 0L) return(invisible(TRUE))

  col_sums <- colSums(m)
  looks_relative <- all(abs(col_sums - 1) < 1e-6) || all(abs(col_sums - 100) < 1e-4)
  msg <- paste0(
    "The feature table holds {length(frac)} non-integer value{?s}",
    if (looks_relative) ", and every sample sums to 1 or 100, so this is a relative abundance table" else "",
    ".\nRarefaction, ANCOM-BC2 and ALDEx2 are defined on counts. Supply counts, ",
    "or pass `expect_counts = FALSE` to state that you know this is not one."
  )
  if (strict) ap_abort(msg) else ap_warn(msg)
  invisible(FALSE)
}

#' @keywords internal
ap_check_no_empty <- function(m, drop_empty_features = TRUE) {
  empty_s <- colnames(m)[colSums(m) == 0]
  ap_assert(
    length(empty_s) == 0L,
    paste0("{length(empty_s)} sample{?s} have zero total counts: ",
           "{paste(utils::head(empty_s, 5), collapse = ', ')}",
           "{if (length(empty_s) > 5) paste0(' and ', length(empty_s) - 5, ' more') else ''}.\n",
           "An empty sample is not a sample. Remove them upstream, where the reason is known.")
  )

  empty_f <- rownames(m)[rowSums(m) == 0]
  if (length(empty_f) > 0L && drop_empty_features) {
    cli::cli_inform(paste0(
      "Dropping {length(empty_f)} feature{?s} with zero counts across all samples. ",
      "They carry no information and inflate the multiple-testing denominator."
    ))
    m <- m[rowSums(m) > 0, , drop = FALSE]
  }
  m
}

# The join. This is the guard that matters most.
#' @keywords internal
ap_check_join <- function(table_ids, meta_ids, max_drop = 0.1,
                          what = c("sample", "feature")) {
  what <- match.arg(what)
  shared <- intersect(table_ids, meta_ids)

  ap_assert(
    length(shared) > 0L,
    paste0("No {what} IDs are shared between the table and the metadata.\n",
           "Table has {length(table_ids)} (e.g. {paste(utils::head(table_ids, 3), collapse = ', ')}).\n",
           "Metadata has {length(meta_ids)} (e.g. {paste(utils::head(meta_ids, 3), collapse = ', ')}).\n",
           "Check that the metadata's first column holds IDs, and that neither side ",
           "was reformatted (numeric IDs read as numbers lose leading zeros).")
  )

  only_table <- setdiff(table_ids, meta_ids)
  only_meta <- setdiff(meta_ids, table_ids)
  drop_frac <- length(only_table) / length(table_ids)

  ap_assert(
    drop_frac <= max_drop,
    paste0("{length(only_table)} of {length(table_ids)} table {what}s ",
           "({round(100 * drop_frac, 1)}%) have no metadata row, above the ",
           "{round(100 * max_drop)}% limit.\n",
           "Examples: {paste(utils::head(only_table, 5), collapse = ', ')}.\n",
           "Raise `max_drop` only once you know why they are missing.")
  )

  if (length(only_table) > 0L) {
    cli::cli_warn(paste0(
      "{length(only_table)} {what}{cli::qty(length(only_table))}{?s} in the table have no metadata row and will be ",
      "dropped: {paste(utils::head(only_table, 5), collapse = ', ')}",
      "{if (length(only_table) > 5) paste0(' and ', length(only_table) - 5, ' more') else ''}."
    ))
  }
  if (length(only_meta) > 0L) {
    cli::cli_inform(paste0(
      "{length(only_meta)} metadata {what}{cli::qty(length(only_meta))}{?s} are absent from the table and will be ",
      "dropped: {paste(utils::head(only_meta, 5), collapse = ', ')}",
      "{if (length(only_meta) > 5) paste0(' and ', length(only_meta) - 5, ' more') else ''}."
    ))
  }

  list(shared = shared, only_table = only_table, only_meta = only_meta)
}

#' @keywords internal
ap_check_tree_covers <- function(tree, feature_ids) {
  ap_assert(inherits(tree, "phylo"), "`tree` must be an `ape::phylo`, not {class(tree)[1]}.")
  missing <- setdiff(feature_ids, tree$tip.label)
  ap_assert(
    length(missing) == 0L,
    paste0("{length(missing)} of {length(feature_ids)} features are absent from the tree: ",
           "{paste(utils::head(missing, 5), collapse = ', ')}.\n",
           "UniFrac and Faith's PD are undefined for a feature with no place on the tree. ",
           "Build the tree from the same representative sequences as the table.")
  )
  invisible(TRUE)
}

#' Re-check an object's internal consistency
#'
#' Asserts that assay dimensions, `colData`, `rowData` and the tree still agree.
#' Called after every operation that changes the sample or feature set, and
#' available to call directly after any manual edit.
#'
#' @param x A `TreeSummarizedExperiment` built by [ap_import()].
#' @return `x`, invisibly.
#' @export
ap_validate <- function(x) {
  counts <- SummarizedExperiment::assay(x, "counts")
  ap_check_count_matrix(counts)

  cd <- SummarizedExperiment::colData(x)
  ap_assert(nrow(cd) == ncol(counts),
            "colData has {nrow(cd)} rows but the assay has {ncol(counts)} samples.")
  ap_assert(identical(rownames(cd), colnames(counts)),
            "colData row names and assay sample IDs disagree in content or order.")

  rd <- SummarizedExperiment::rowData(x)
  ap_assert(nrow(rd) == nrow(counts),
            "rowData has {nrow(rd)} rows but the assay has {nrow(counts)} features.")
  ap_assert(identical(rownames(rd), rownames(counts)),
            "rowData row names and assay feature IDs disagree in content or order.")

  tree <- tryCatch(TreeSummarizedExperiment::rowTree(x), error = function(e) NULL)
  if (!is.null(tree)) ap_check_tree_covers(tree, rownames(counts))

  invisible(x)
}
