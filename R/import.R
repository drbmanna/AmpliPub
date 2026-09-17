#' Import an amplicon dataset
#'
#' The entry point. Takes a feature table and sample metadata, optionally a
#' phylogeny and taxonomy, runs the join and consistency guards, and returns a
#' `TreeSummarizedExperiment` carrying an AmpliPub provenance log.
#'
#' AmpliPub is pipeline-agnostic. Anything that yields a feature table works:
#' QIIME 2 `.qza`, BIOM (v1 JSON or v2 HDF5), a plain TSV, or an in-memory
#' matrix. QIIME 2 artifacts are read directly, without QIIME 2 installed.
#'
#' The returned object is a plain `TreeSummarizedExperiment`, not a wrapper
#' class, so `mia`, `ANCOMBC` and the rest of the Bioconductor stack accept it
#' as-is. The provenance log lives in its `metadata()`; see [ap_summary()].
#'
#' @section Sample IDs are read as text:
#' Metadata is read with every column as character, then types are inferred
#' column by column, because sample IDs are the one column that must never be
#' guessed. Baxter 2016 IDs look like `2003650`. Read as numbers they survive;
#' an ID like `007` does not. A silently reformatted ID does not fail the join
#' loudly, it fails it quietly for the subset that happened to look numeric.
#'
#' @param table Feature table: path to a `.qza`, `.biom` or `.tsv`, or a numeric
#'   matrix with features as rows.
#' @param metadata Sample metadata: path to a TSV/CSV whose first column holds
#'   sample IDs, or a data frame with sample IDs in the row names or first column.
#' @param tree Optional phylogeny: path to a `.qza` or Newick file, or an
#'   `ape::phylo`. Required for UniFrac and Faith's PD.
#' @param taxonomy Optional taxonomy: path to a `.qza` or TSV, or a data frame
#'   with `feature_id` and `taxon` columns.
#' @param features_are_rows Orientation of a TSV or matrix input. Default `TRUE`.
#' @param expect_counts Warn when the table is not integer counts. Default `TRUE`.
#' @param max_drop Largest fraction of table samples allowed to lack a metadata
#'   row before import fails. Default `0.1`.
#' @param drop_empty_features Remove features with zero counts everywhere, which
#'   is normal after subsetting. Default `TRUE`.
#'
#' @return A `TreeSummarizedExperiment` with a `counts` assay.
#' @export
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(200, 20), nrow = 20,
#'                  dimnames = list(paste0("ASV", 1:20), paste0("S", 1:10)))
#' meta <- data.frame(group = rep(c("a", "b"), each = 5),
#'                    row.names = paste0("S", 1:10))
#' x <- ap_import(counts, meta)
#' ap_summary(x)
ap_import <- function(table,
                      metadata,
                      tree = NULL,
                      taxonomy = NULL,
                      features_are_rows = TRUE,
                      expect_counts = TRUE,
                      max_drop = 0.1,
                      drop_empty_features = TRUE) {

  log <- ap_new_log()

  counts <- ap_load_table(table, features_are_rows = features_are_rows)
  log <- ap_log_input(log, "table", table,
                      list(features = nrow(counts), samples = ncol(counts)))

  ap_check_count_matrix(counts)
  if (expect_counts) ap_check_counts_are_integers(counts, strict = FALSE)

  meta <- ap_load_metadata(metadata)
  log <- ap_log_input(log, "metadata", metadata,
                      list(samples = nrow(meta), variables = ncol(meta)))

  join <- ap_check_join(colnames(counts), rownames(meta), max_drop = max_drop,
                        what = "sample")
  # Order by the table, not by the metadata, so the assay is never permuted
  # relative to what the user handed in.
  keep <- colnames(counts)[colnames(counts) %in% join$shared]
  counts <- counts[, keep, drop = FALSE]
  meta <- meta[keep, , drop = FALSE]

  counts <- ap_check_no_empty(counts, drop_empty_features = drop_empty_features)

  # Features go into a fixed order, because upstream order is not stable and some results
  # depend on it. The same Snakemake config run twice produced the same 829 features in a
  # different order, and ALDEx2 draws its Monte Carlo instances feature by feature from the
  # RNG stream, so an identical seed on a differently ordered table gave different draws per
  # feature. A seed cannot fix that; the order has to be imposed. `method = "radix"` sorts
  # in C collation rather than the session locale, so the order is the same on every
  # machine. Samples are deliberately left in the table's own order, per the join above.
  counts <- counts[order(rownames(counts), method = "radix"), , drop = FALSE]

  row_data <- S4Vectors::DataFrame(row.names = rownames(counts))
  if (!is.null(taxonomy)) {
    tax <- ap_load_taxonomy(taxonomy)
    log <- ap_log_input(log, "taxonomy", taxonomy, list(features = nrow(tax)))
    row_data <- ap_attach_taxonomy(rownames(counts), tax)
  }

  phy <- NULL
  if (!is.null(tree)) {
    phy <- ap_load_tree(tree)
    log <- ap_log_input(log, "tree", tree, list(tips = length(phy$tip.label)))
    ap_check_tree_covers(phy, rownames(counts))
    # Tips not in the table are dropped so downstream UniFrac is computed on the
    # tree the table actually spans.
    extra <- setdiff(phy$tip.label, rownames(counts))
    if (length(extra) > 0L) {
      cli::cli_inform("Pruning {length(extra)} tree tip{?s} absent from the table.")
      phy <- ape::drop.tip(phy, extra)
    }
  }

  x <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(meta, row.names = rownames(meta)),
    rowData = row_data,
    rowTree = phy
  )

  x <- ap_set_log(x, log)
  x <- ap_log_step(x, "import", list(
    features = nrow(counts),
    samples = ncol(counts),
    total_counts = sum(counts),
    min_depth = min(colSums(counts)),
    max_depth = max(colSums(counts)),
    dropped_no_metadata = length(join$only_table),
    metadata_without_sample = length(join$only_meta)
  ))
  ap_validate(x)
  x
}

#' @keywords internal
ap_load_table <- function(table, features_are_rows = TRUE) {
  if (is.matrix(table) || is.data.frame(table)) {
    m <- as.matrix(table)
    storage.mode(m) <- "double"
  } else {
    ap_assert(is.character(table) && length(table) == 1L,
              "`table` must be a file path or a matrix, not {class(table)[1]}.")
    ap_assert(file.exists(table), "Feature table not found: {table}")
    ext <- tolower(tools::file_ext(table))
    m <- switch(
      ext,
      qza = {
        obj <- ap_read_qza(table)
        ap_assert(
          startsWith(attr(obj, "qza_type") %||% "", "FeatureTable"),
          paste0("{table} is a QIIME 2 artifact of type `{attr(obj, 'qza_type')}`, ",
                 "not a FeatureTable.")
        )
        obj
      },
      biom = ap_read_biom_file(table),
      ap_read_table_tsv(table)
    )
  }
  if (!features_are_rows) m <- t(m)
  m
}

#' @keywords internal
ap_read_table_tsv <- function(path) {
  lines <- readLines(path, n = 2L, warn = FALSE)
  # `biom convert --to-tsv` writes a `# Constructed from biom file` banner above
  # the header.
  skip <- if (length(lines) > 0L && startsWith(lines[1], "# ")) 1L else 0L
  sep <- if (tolower(tools::file_ext(path)) == "csv") "," else "\t"
  df <- utils::read.delim(path, sep = sep, header = TRUE, row.names = 1,
                          skip = skip, quote = "", check.names = FALSE,
                          comment.char = "")
  m <- as.matrix(df)
  ap_assert(is.numeric(m),
            paste0("Feature table {path} has non-numeric columns after the ID column. ",
                   "If taxonomy is stored as a trailing column, remove it and pass it ",
                   "through the `taxonomy` argument."))
  storage.mode(m) <- "double"
  m
}

#' @keywords internal
ap_load_metadata <- function(metadata) {
  if (is.data.frame(metadata)) {
    df <- metadata
    if (is.null(rownames(df)) || identical(rownames(df), as.character(seq_len(nrow(df))))) {
      ap_assert(ncol(df) >= 1L, "Metadata data frame has no columns.")
      rownames(df) <- as.character(df[[1]])
      df <- df[, -1, drop = FALSE]
    }
    rownames(df) <- as.character(rownames(df))
    return(df)
  }

  ap_assert(is.character(metadata) && length(metadata) == 1L,
            "`metadata` must be a file path or a data frame, not {class(metadata)[1]}.")
  ap_assert(file.exists(metadata), "Metadata file not found: {metadata}")

  sep <- if (tolower(tools::file_ext(metadata)) == "csv") "," else "\t"
  raw <- utils::read.delim(metadata, sep = sep, header = TRUE, quote = "",
                           check.names = FALSE, comment.char = "",
                           colClasses = "character", na.strings = c("NA"))
  ap_assert(ncol(raw) >= 2L,
            "Metadata file {metadata} has {ncol(raw)} column{?s}; expected an ID column plus at least one variable.")

  # QIIME 2 metadata may carry a directive row (`#q2:types`) under the header.
  if (grepl("^#", raw[[1]][1])) raw <- raw[-1, , drop = FALSE]

  ids <- trimws(raw[[1]])
  ap_assert(!anyDuplicated(ids),
            "Metadata file {metadata} has duplicate sample IDs: {paste(utils::head(ids[duplicated(ids)], 5), collapse = ', ')}.")
  ap_assert(all(nzchar(ids)), "Metadata file {metadata} has blank sample IDs.")

  df <- raw[, -1, drop = FALSE]
  rownames(df) <- ids
  ap_infer_types(df)
}

# Types are inferred per column AFTER the IDs are safely held as text. A column
# is numeric only if every non-missing value parses as a number; anything else
# stays character.
#' @keywords internal
ap_infer_types <- function(df) {
  for (j in seq_along(df)) {
    v <- trimws(df[[j]])
    v[!nzchar(v)] <- NA_character_
    num <- suppressWarnings(as.numeric(v))
    if (all(is.na(num) == is.na(v)) && any(!is.na(num))) {
      df[[j]] <- num
    } else {
      df[[j]] <- v
    }
  }
  df
}

#' @keywords internal
ap_load_tree <- function(tree) {
  if (inherits(tree, "phylo")) return(tree)
  ap_assert(is.character(tree) && length(tree) == 1L,
            "`tree` must be a file path or an `ape::phylo`, not {class(tree)[1]}.")
  ap_assert(file.exists(tree), "Tree file not found: {tree}")
  if (tolower(tools::file_ext(tree)) == "qza") {
    obj <- ap_read_qza(tree)
    ap_assert(inherits(obj, "phylo"),
              "{tree} is a QIIME 2 artifact of type `{attr(obj, 'qza_type')}`, not a Phylogeny.")
    return(obj)
  }
  ape::read.tree(tree)
}

#' @keywords internal
ap_load_taxonomy <- function(taxonomy) {
  if (is.data.frame(taxonomy)) {
    ap_assert(all(c("feature_id", "taxon") %in% names(taxonomy)),
              "A taxonomy data frame needs `feature_id` and `taxon` columns.")
    return(taxonomy)
  }
  ap_assert(is.character(taxonomy) && length(taxonomy) == 1L,
            "`taxonomy` must be a file path or a data frame, not {class(taxonomy)[1]}.")
  ap_assert(file.exists(taxonomy), "Taxonomy file not found: {taxonomy}")
  if (tolower(tools::file_ext(taxonomy)) == "qza") {
    obj <- ap_read_qza(taxonomy)
    ap_assert(is.data.frame(obj),
              "{taxonomy} is a QIIME 2 artifact of type `{attr(obj, 'qza_type')}`, not FeatureData[Taxonomy].")
    return(obj)
  }
  ap_read_q2_taxonomy(taxonomy)
}

#' @keywords internal
ap_attach_taxonomy <- function(feature_ids, tax) {
  idx <- match(feature_ids, tax$feature_id)
  n_missing <- sum(is.na(idx))
  if (n_missing > 0L) {
    ap_warn(paste0(
      "{n_missing} of {length(feature_ids)} features have no taxonomy assignment and ",
      "will read as Unassigned: {paste(utils::head(feature_ids[is.na(idx)], 5), collapse = ', ')}."
    ))
  }
  taxon <- tax$taxon[idx]
  ranks <- ap_parse_taxonomy(ifelse(is.na(taxon), "", taxon))
  out <- S4Vectors::DataFrame(taxon = taxon, ranks, row.names = feature_ids)
  if ("confidence" %in% names(tax)) out$confidence <- tax$confidence[idx]
  out
}
