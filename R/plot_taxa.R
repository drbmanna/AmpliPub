#' Stacked composition bar chart
#'
#' Relative abundance of the most abundant taxa, one bar per sample or one per
#' group.
#'
#' The caption states how many taxa were pooled into `Other` and how much
#' abundance they carry. A stacked bar chart showing 15 genera out of 250 looks
#' like a complete picture, and whether it is depends entirely on how large
#' `Other` is.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()], or a data frame
#'   from [ap_top_taxa()].
#' @param rank Rank to collapse to. Default `"genus"`.
#' @param n Number of taxa to show. Default `15`.
#' @param group Optional metadata variable to facet by.
#' @param mode `"sample"` (default) draws one bar per sample; `"group"` draws
#'   the mean composition per group.
#' @param order_by Taxon whose abundance orders the samples within each facet.
#'   Defaults to the most abundant. Ignored when `mode = "group"`.
#'
#' @return A ggplot object.
#' @export
ap_plot_taxa_bar <- function(x, rank = "genus", n = 15L, group = NULL,
                             mode = c("sample", "group"), order_by = NULL) {
  mode <- match.arg(mode)
  df <- if (is.data.frame(x)) x else ap_top_taxa(x, n = n, group = group, rank = rank)
  group <- group %||% attr(df, "group")
  # aggregate() drops attributes, and the pooling counts live there, so the
  # caption is built before any reshaping.
  caption <- ap_taxa_caption(df, rank)

  if (mode == "group") {
    ap_assert(!is.null(group) && "group" %in% names(df),
              "`mode = \"group\"` needs a `group` variable to average within.")
    taxon_levels <- levels(df$taxon)
    df <- stats::aggregate(relative_abundance ~ taxon + group, data = df, FUN = mean)
    df$taxon <- factor(as.character(df$taxon), levels = taxon_levels)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$relative_abundance,
                                          fill = .data$taxon)) +
      ggplot2::geom_col(width = 0.75, colour = "white", linewidth = 0.15) +
      ggplot2::labs(x = group)
  } else {
    ref <- order_by %||% levels(df$taxon)[1]
    ord <- df[df$taxon == ref, c("sample_id", "relative_abundance")]
    ord <- ord[order(-ord$relative_abundance), ]
    df$sample_id <- factor(df$sample_id, levels = ord$sample_id)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sample_id,
                                          y = .data$relative_abundance,
                                          fill = .data$taxon)) +
      ggplot2::geom_col(width = 1) +
      ggplot2::labs(x = paste0("Sample (ordered by ", ref, ")")) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    if (!is.null(group) && "group" %in% names(df)) {
      p <- p + ggplot2::facet_wrap(~ group, scales = "free_x")
    }
  }

  p +
    ap_scale_fill_taxa(levels(df$taxon), name = rank) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(y = "Relative abundance", caption = caption) +
    ap_theme(grid = "none") +
    ggplot2::theme(axis.text.x = if (mode == "sample") ggplot2::element_blank()
                   else ggplot2::element_text(),
                   axis.ticks.x = if (mode == "sample") ggplot2::element_blank()
                   else ggplot2::element_line())
}

#' @keywords internal
ap_taxa_caption <- function(df, rank) {
  pooled <- attr(df, "n_pooled")
  if (is.null(pooled) || pooled == 0L) {
    return(paste0("All ", rank, "-level taxa shown."))
  }
  ap_wrap(sprintf(
    paste0("%d %s-level taxa shown; %d further taxa pooled as Other, carrying %.1f%% ",
           "of mean relative abundance. Taxa unassigned at %s keep the label of the ",
           "deepest rank they reached."),
    attr(df, "n_kept"), rank, pooled,
    100 * attr(df, "pooled_mean_abundance"), rank))
}

# ggplot2 does not wrap caption text, so a long caption runs off the device.
#' @keywords internal
ap_wrap <- function(text, width = 95L) {
  paste(strwrap(text, width = width), collapse = "\n")
}

#' Heatmap of taxon abundance
#'
#' Log-scaled relative abundance, taxa by samples, optionally split by a
#' metadata variable.
#'
#' Log scale because abundance spans orders of magnitude and a linear colour
#' scale renders everything below the top few taxa as the same shade of nothing.
#' The pseudocount that makes the log defined is stated in the caption.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param rank Rank to collapse to. Default `"genus"`.
#' @param n Number of taxa to show. Default `15`.
#' @param group Optional metadata variable to facet by.
#' @param pseudocount Added before the log. Default `1e-5`.
#'
#' @return A ggplot object.
#' @export
ap_plot_taxa_heatmap <- function(x, rank = "genus", n = 15L, group = NULL,
                                 pseudocount = 1e-5) {
  df <- ap_top_taxa(x, n = n, group = group, rank = rank)
  df <- df[df$taxon != "Other", , drop = FALSE]
  df$taxon <- droplevels(df$taxon)
  df$log_abundance <- log10(df$relative_abundance + pseudocount)

  if (!is.null(group) && "group" %in% names(df)) {
    ord <- stats::aggregate(relative_abundance ~ sample_id, data = df, FUN = sum)
    df$sample_id <- factor(df$sample_id, levels = ord$sample_id[order(-ord$relative_abundance)])
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$sample_id,
                                        y = .data$taxon,
                                        fill = .data$log_abundance)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(
      name = expression(log[10] * " rel. abundance"), option = "mako", direction = -1
    ) +
    ggplot2::scale_y_discrete(limits = rev(levels(df$taxon))) +
    ggplot2::labs(x = "Sample", y = NULL,
                  caption = sprintf("Pseudocount %g added before the log.", pseudocount)) +
    ap_theme(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())

  if (!is.null(group) && "group" %in% names(df)) {
    p <- p + ggplot2::facet_wrap(~ group, scales = "free_x", nrow = 1)
  }
  p
}
