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
#' @param publication Draw the publication version: mean composition per group,
#'   rank prefixes removed, taxon names in italics, no caption. Needs `group`.
#'   The caption's content goes to [ap_taxa_legend()] instead.
#' @param pub [ap_pub_options()] for `publication = TRUE`.
#'
#' @return A ggplot object.
#' @export
ap_plot_taxa_bar <- function(x, rank = "genus", n = 15L, group = NULL,
                             mode = c("sample", "group"), order_by = NULL,
                             publication = FALSE, pub = ap_pub_options()) {
  mode <- match.arg(mode)
  df <- if (is.data.frame(x)) x else ap_top_taxa(x, n = n, group = group, rank = rank)
  group <- group %||% attr(df, "group")
  if (publication) return(ap_plot_taxa_bar_pub(df, group, pub, rank))
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
#' @param publication Draw the publication version, which compares groups rather than
#'   showing samples: one column per group, each cell the group's mean CLR for that
#'   taxon minus a baseline, on a diverging scale centred on 0. Needs `group`. See the
#'   section below.
#' @param reference For `publication = TRUE`: a level of `group` to use as the baseline.
#'   `NULL` (default) uses the average of the group means instead, and every group gets
#'   a column. When a reference is given its column would be all zeros, so it is not
#'   drawn; the colour bar title and [ap_taxa_legend()] name the reference instead.
#' @param clr_pseudocount For `publication = TRUE`: the value that replaces zeros before
#'   the CLR, passed to [ap_normalize()]. Default `0.5`.
#' @param pub [ap_pub_options()] for `publication = TRUE`.
#'
#' @section The publication heatmap:
#' Relative abundances are compositional: one taxon growing shrinks every other taxon's
#' share, so a difference in proportions between groups can be produced by other taxa.
#' The publication heatmap therefore works on the centred log-ratio (CLR) scale, computed
#' on every taxon at `rank`, not only the ones shown. Each cell is the group's mean CLR
#' minus a baseline: the reference group's mean when `reference` is given (its column
#' is then all zeros and is dropped), otherwise the unweighted average of the group
#' means, so no group has to be the reference and a large group does not pull the
#' centre toward itself. The row is centred but not scaled: a z-score would stretch
#' every row to the same range, however small the real difference, and with few groups
#' that makes noise look like a finding. The figure is descriptive; differential
#' abundance is tested by [ap_da()].
#'
#' @return A ggplot object.
#' @export
ap_plot_taxa_heatmap <- function(x, rank = "genus", n = 15L, group = NULL,
                                 pseudocount = 1e-5, publication = FALSE,
                                 reference = NULL, clr_pseudocount = 0.5,
                                 pub = ap_pub_options()) {
  df <- ap_top_taxa(x, n = n, group = group, rank = rank)
  if (publication) {
    return(ap_plot_taxa_heatmap_pub(x, df, group, rank, reference, clr_pseudocount, pub))
  }
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

# Publication versions -----------------------------------------------------------------------

# Taxon names as the reader sees them: the rank prefix removed and every suffix kept
# (GTDB letter suffixes and Greengenes2 numeric identifiers are part of the name, and
# removing them would merge distinct lineages), in italics. A label that fell back to a
# higher rank keeps its "(unassigned genus)" note upright. "Other" is not a taxon and is
# not italic. Returns plotmath expressions named by the original label.
#' @keywords internal
ap_pub_taxon_labels <- function(taxa) {
  taxa <- as.character(taxa)
  shown <- sub("^[a-z]__", "", taxa)
  ap_assert(!anyDuplicated(shown),
            "Removing rank prefixes would give two taxa the same label: {paste(unique(shown[duplicated(shown)]), collapse = ', ')}.")
  fallback <- regmatches(shown, regexec("^(.*) \\((unassigned [a-z]+)\\)$", shown))
  out <- lapply(seq_along(shown), function(i) {
    s <- shown[i]
    if (s == "Other") return("Other")
    m <- fallback[[i]]
    if (length(m) == 3L) {
      note <- paste0("(", m[3], ")")
      bquote(italic(.(m[2])) ~ .(note))
    } else {
      bquote(italic(.(s)))
    }
  })
  stats::setNames(as.expression(out), taxa)
}

# More colours than any ggsci journal palette holds: stacked bars routinely show 15 to 20
# taxa. d3 category20 without its two greys, which would be confused with Other, then
# IGV colours for anything beyond. Not chosen for colour vision deficiency.
# "genus" -> "Genus", the legend or axis title naming what the taxa are.
#' @keywords internal
ap_rank_title <- function(rank) paste0(toupper(substr(rank, 1, 1)), substring(rank, 2))

#' @keywords internal
ap_pub_taxa_palette <- function(n) {
  ap_pub_require("ggsci")
  is_grey <- function(v) {
    rgb <- grDevices::col2rgb(v)
    apply(rgb, 2, function(z) diff(range(z)) < 16)
  }
  pal <- substr(ggsci::pal_d3("category20")(20), 1, 7)
  pal <- pal[!is_grey(pal)]
  if (n > length(pal)) {
    extra <- substr(ggsci::pal_igv("default")(51), 1, 7)
    extra <- extra[!is_grey(extra) & !toupper(extra) %in% toupper(pal)]
    pal <- c(pal, extra)
  }
  ap_assert(n <= length(pal),
            "{n} taxa need more colours than are available ({length(pal)}). Lower `n`.")
  pal[seq_len(n)]
}

#' @keywords internal
ap_plot_taxa_bar_pub <- function(df, group, pub, rank) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  ap_assert(!is.null(group) && "group" %in% names(df),
            "The publication composition figure shows group means, so it needs `group`.")
  taxon_levels <- levels(df$taxon)
  df <- stats::aggregate(relative_abundance ~ taxon + group, data = df, FUN = mean)
  df$taxon <- factor(as.character(df$taxon), levels = taxon_levels)
  df$group <- ap_pub_levels(df$group, pub)

  real <- setdiff(taxon_levels, "Other")
  values <- stats::setNames(ap_pub_taxa_palette(length(real)), real)
  if ("Other" %in% taxon_levels) values <- c(values, Other = "#D9D9D9")
  labels <- ap_pub_taxon_labels(taxon_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$relative_abundance,
                                        fill = .data$taxon)) +
    ggplot2::geom_col(width = 0.7, colour = "white", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = values, breaks = taxon_levels,
                               labels = labels[taxon_levels], name = ap_rank_title(rank)) +
    ggplot2::scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(x = ap_pub_label(group, pub), y = "Mean relative abundance") +
    ap_theme_pub(legend = "right") +
    ggplot2::theme(legend.key.size = ggplot2::unit(3, "mm"),
                   legend.text = ggplot2::element_text(size = 6),
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
  attr(p, "ap_pub_size") <- list(width = "onehalf",
                                 height = max(60, 3.3 * length(taxon_levels) + 12))
  p
}

# Group mean CLR minus the baseline, for the taxa `df` selected. Returns taxa by groups.
#' @keywords internal
ap_taxa_clr_contrast <- function(x, df, group, rank, reference, clr_pseudocount) {
  collapsed <- ap_collapse(x, rank = rank)
  clr <- ap_normalize(collapsed, method = "clr", pseudocount = clr_pseudocount)
  g <- SummarizedExperiment::colData(collapsed)[[group]][
    match(colnames(clr), colnames(collapsed))]
  keep <- !is.na(g)
  g <- factor(g[keep])
  shown <- setdiff(levels(df$taxon), "Other")
  ap_assert(all(shown %in% rownames(clr)),
            "Taxa in the figure are missing from the collapsed table: {paste(setdiff(shown, rownames(clr)), collapse = ', ')}.")
  means <- vapply(levels(g), function(lv) rowMeans(clr[shown, keep, drop = FALSE][, g == lv, drop = FALSE]),
                  numeric(length(shown)))
  means <- matrix(means, nrow = length(shown), dimnames = list(shown, levels(g)))
  if (is.null(reference)) {
    means - rowMeans(means)
  } else {
    ap_assert(reference %in% levels(g),
              "`reference` must be a level of `{group}`: {paste(levels(g), collapse = ', ')}.")
    out <- means - means[, reference]
    out[, colnames(out) != reference, drop = FALSE]
  }
}

#' @keywords internal
ap_plot_taxa_heatmap_pub <- function(x, df, group, rank, reference, clr_pseudocount, pub) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  ap_assert(!is.null(group) && "group" %in% names(df),
            "The publication heatmap compares groups, so it needs `group`.")
  m <- ap_taxa_clr_contrast(x, df, group, rank, reference, clr_pseudocount)
  taxon_levels <- rownames(m)
  long <- data.frame(
    taxon = factor(rep(taxon_levels, times = ncol(m)), levels = taxon_levels),
    group = rep(colnames(m), each = nrow(m)),
    value = as.vector(m)
  )
  long$group <- ap_pub_levels(factor(long$group, levels = colnames(m)), pub)
  lim <- max(abs(long$value))
  labels <- ap_pub_taxon_labels(taxon_levels)
  fill_title <- if (is.null(reference)) "Mean CLR relative to the group average" else
    paste0("Mean CLR relative to ", levels(ap_pub_levels(factor(reference), pub)))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$group, y = .data$taxon,
                                          fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(name = fill_title, low = "#2166AC", mid = "white",
                                  high = "#B2182B", midpoint = 0, limits = c(-lim, lim)) +
    ggplot2::scale_y_discrete(limits = rev(taxon_levels), labels = labels[rev(taxon_levels)]) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = ap_pub_label(group, pub), y = ap_rank_title(rank)) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    ap_theme_pub(legend = "bottom") +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size = 6),
                   legend.key.height = ggplot2::unit(2.2, "mm"),
                   legend.key.width = ggplot2::unit(6, "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
  attr(p, "ap_pub_size") <- list(width = "single",
                                 height = max(55, 3.4 * length(taxon_levels) + 30))
  p
}

#' Legend text for a publication composition figure
#'
#' What the report composition figures write in their caption, for the figure
#' legend: how the taxa were chosen and how many are shown, what `Other` carries,
#' how unassigned features are labelled, and the number of samples per group.
#' Written next to the figure by [ap_save_figure()].
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param rank Rank the figure was drawn at.
#' @param group Metadata variable the figure groups by.
#' @param n Taxa per group, as passed to the plot. Default `15`.
#' @param type `"bar"` or `"heatmap"`.
#' @param reference For `type = "heatmap"`, the reference level, or `NULL` when the
#'   heatmap is centred on the average of the group means.
#' @param clr_pseudocount For `type = "heatmap"`, the zero replacement used before the CLR.
#' @param pub [ap_pub_options()], for the group names (`level_labels`).
#' @return A character vector, one line per entry.
#' @export
ap_taxa_legend <- function(x, rank, group, n = 15L, type = c("bar", "heatmap"),
                           reference = NULL, clr_pseudocount = 0.5, pub = ap_pub_options()) {
  type <- match.arg(type)
  df <- ap_top_taxa(x, n = n, group = group, rank = rank)
  sizes <- table(df$group[!duplicated(df$sample_id)])
  pooled <- attr(df, "n_pooled")
  c(
    if (type == "bar") {
      "Mean relative abundance per group; counts converted to proportions within each sample."
    } else {
      paste0(sprintf(paste0("Mean centred log-ratio (CLR) per group, relative to %s. CLR computed ",
                            "on all %s-level taxa after replacing zeros with %g. "),
                     if (is.null(reference)) "the average of the group means"
                     else paste0("the ", ap_pub_level_text(reference, pub),
                                 " group (the reference, not drawn as a column)"),
                     rank, clr_pseudocount),
             "Descriptive; differential abundance is tested separately.")
    },
    sprintf(paste0("Taxa shown: the %d most abundant %s-level taxa by mean relative abundance ",
                   "within each %s group, combined across groups (%d taxa)."),
            as.integer(n), rank, group, attr(df, "n_kept")),
    if (is.null(pooled) || pooled == 0L) {
      "No further taxa."
    } else {
      sprintf("%d further taxa carry %.1f%% of mean relative abundance%s.", pooled,
              100 * attr(df, "pooled_mean_abundance"),
              if (type == "bar") " and are pooled as Other" else " and are not shown")
    },
    sprintf("Features unassigned at %s are labelled with the deepest rank they reached.", rank),
    paste0("n = ", paste(sprintf("%d %s", as.integer(sizes), ap_pub_level_text(names(sizes), pub)),
                         collapse = ", "), ".")
  )
}
