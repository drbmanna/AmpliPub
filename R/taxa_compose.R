# Composition: collapsing to a rank and picking the taxa worth showing.
#
# Two things here are usually done silently and both change what the figure
# says. Collapsing to a rank has to decide what to do with features unassigned
# at that rank, and a top-N bar chart has to decide what "top" means and where
# everything else goes. Dropping the unassigned features quietly inflates every
# remaining proportion, because the denominator shrinks.

#' Collapse a table to a taxonomic rank
#'
#' Sums counts across features sharing a label at `rank`. Features unassigned
#' at that rank fall back up the hierarchy and keep a label that says so (see
#' [ap_taxon_label()]), rather than being dropped.
#'
#' Dropping them is the common alternative and it is not neutral. On this
#' dataset 91 of 654 features have no genus, carrying real abundance; removing
#' them would raise every reported genus proportion without saying so.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()] with taxonomy.
#' @param rank Rank to collapse to.
#' @param drop_unassigned Remove features with no assignment at any rank.
#'   Default `FALSE`.
#'
#' @return A `TreeSummarizedExperiment` with one row per taxon. The tree is
#'   dropped, since a collapsed taxon has no single place on it.
#' @export
ap_collapse <- function(x, rank = "genus", drop_unassigned = FALSE) {
  rd <- as.data.frame(SummarizedExperiment::rowData(x))
  ap_assert(rank %in% colnames(rd),
            paste0("Rank `{rank}` is not present. This object was imported without ",
                   "taxonomy, or at a shallower depth. Available: ",
                   "{paste(intersect(ap_rank_names(), colnames(rd)), collapse = ', ')}."))

  counts <- SummarizedExperiment::assay(x, "counts")
  label <- ap_taxon_label(rd, rank)

  if (drop_unassigned) {
    keep <- label != "Unassigned"
    cli::cli_inform(paste0(
      "Dropping {sum(!keep)} unassigned feature{?s} carrying ",
      "{sprintf('%.2f%%', 100 * sum(counts[!keep, ]) / sum(counts))} of all reads. ",
      "Every remaining proportion is now relative to the smaller total."
    ))
    counts <- counts[keep, , drop = FALSE]
    label <- label[keep]
  }

  collapsed <- rowsum(counts, group = label, reorder = TRUE)
  fallback <- grepl("(unassigned ", label, fixed = TRUE)

  out <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = collapsed),
    colData = SummarizedExperiment::colData(x),
    rowData = S4Vectors::DataFrame(
      taxon = rownames(collapsed),
      rank = rank,
      n_features = as.vector(table(label)[rownames(collapsed)]),
      row.names = rownames(collapsed)
    )
  )
  out <- ap_set_log(out, ap_get_log(x))
  ap_log_step(out, "collapse", list(
    rank = rank,
    features_in = nrow(counts),
    taxa_out = nrow(collapsed),
    unassigned_at_rank = sum(fallback) + sum(label == "Unassigned"),
    reads_unchanged = sum(collapsed) == sum(counts)
  ))
}

#' Most abundant taxa
#'
#' Ranks taxa by mean relative abundance across samples and pools the rest as
#' `Other`.
#'
#' Mean relative abundance, not total counts. Ranking on totals lets a handful
#' of deeply sequenced samples decide which taxa appear, which is a property of
#' the library preparation rather than of the biology.
#'
#' @param x A `TreeSummarizedExperiment`, usually from [ap_collapse()].
#' @param n Number of taxa to keep. Default `15`.
#' @param group Optional metadata variable. When given, a taxon qualifies if it
#'   makes the top `n` in any group, so a taxon abundant in a small group is
#'   not hidden by a large one.
#' @param rank Rank to collapse to first, when `x` has not been collapsed.
#'
#' @return A data frame: `sample_id`, `taxon`, `relative_abundance`, and `group`
#'   when given. `taxon` is a factor ordered by overall abundance with `Other`
#'   last.
#' @export
ap_top_taxa <- function(x, n = 15L, group = NULL, rank = NULL) {
  if (!is.null(rank)) x <- ap_collapse(x, rank)
  counts <- SummarizedExperiment::assay(x, "counts")
  rel <- sweep(counts, 2, colSums(counts), "/")
  meta <- as.data.frame(SummarizedExperiment::colData(x))

  overall <- sort(rowMeans(rel), decreasing = TRUE)
  keep <- names(overall)[seq_len(min(n, length(overall)))]

  if (!is.null(group)) {
    ap_assert(group %in% names(meta), "Variable `{group}` is not in the sample metadata.")
    g <- meta[[group]][match(colnames(rel), rownames(meta))]
    for (lv in unique(g[!is.na(g)])) {
      cols <- which(!is.na(g) & g == lv)
      if (length(cols) == 0L) next
      means <- sort(rowMeans(rel[, cols, drop = FALSE]), decreasing = TRUE)
      keep <- union(keep, names(means)[seq_len(min(n, length(means)))])
    }
  }

  pooled <- setdiff(rownames(rel), keep)
  long <- data.frame(
    sample_id = rep(colnames(rel), each = length(keep)),
    taxon = rep(keep, times = ncol(rel)),
    relative_abundance = as.vector(rel[keep, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  if (length(pooled) > 0L) {
    other <- colSums(rel[pooled, , drop = FALSE])
    long <- rbind(long, data.frame(
      sample_id = colnames(rel), taxon = "Other",
      relative_abundance = as.vector(other), stringsAsFactors = FALSE
    ))
  }

  levels_ordered <- c(names(overall)[names(overall) %in% keep], "Other")
  levels_ordered <- levels_ordered[levels_ordered %in% unique(long$taxon)]
  long$taxon <- factor(long$taxon, levels = levels_ordered)

  if (!is.null(group)) {
    long$group <- meta[[group]][match(long$sample_id, rownames(meta))]
  }

  attr(long, "n_kept") <- length(keep)
  attr(long, "n_pooled") <- length(pooled)
  attr(long, "pooled_mean_abundance") <- if (length(pooled) > 0L) sum(overall[pooled]) else 0
  attr(long, "group") <- group
  long
}
