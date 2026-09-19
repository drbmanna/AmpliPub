#' Plot an ordination
#'
#' PCoA or NMDS with the numbers that make it readable: variance explained on
#' the axes, stress in the corner for NMDS, and the PERMANOVA result with its
#' dispersion verdict written onto the panel.
#'
#' The dispersion verdict is on the figure rather than in the caption on
#' purpose. An ellipse plot with `p = 0.001` reads as "these groups are
#' different"; if the dispersion test also fired, that reading is wrong, and the
#' figure is where the reader will form it.
#'
#' @param ord An `ap_ordination` from [ap_ordinate()].
#' @param group Metadata variable to colour by.
#' @param shape Optional second metadata variable mapped to point shape.
#' @param permanova Optional `ap_permanova` result to annotate with.
#' @param ellipse Draw a 95% confidence ellipse per group. Default `TRUE`.
#' @param axes Which two axes to plot. Default `c(1, 2)`.
#' @param point_size Point size. Default `1.6`.
#' @param publication Draw the publication version: no subtitle or caption,
#'   [ap_theme_pub()], ggsci colours, ellipses filled at low opacity, and the
#'   PERMANOVA and dispersion numbers on the panel without the verdict. The
#'   verdict, its explanation and the projection caveat go to the legend text
#'   from [ap_ordination_legend()]. Default `FALSE`.
#' @param pub [ap_pub_options()] for `publication = TRUE`.
#'
#' @return A ggplot object.
#' @export
ap_plot_ordination <- function(ord, group = NULL, shape = NULL, permanova = NULL,
                               ellipse = TRUE, axes = c(1, 2), point_size = 1.6,
                               publication = FALSE, pub = ap_pub_options()) {
  ap_assert(inherits(ord, "ap_ordination"),
            "`ord` must come from `ap_ordinate()`, not {class(ord)[1]}.")
  ap_assert(length(axes) == 2L && all(axes <= ncol(ord$coords)),
            "`axes` must name two axes that exist; this ordination has {ncol(ord$coords)}.")

  df <- as.data.frame(ord$coords[, axes, drop = FALSE])
  names(df) <- c("x", "y")
  df$sample_id <- rownames(ord$coords)

  meta <- ord$metadata
  if (!is.null(group)) {
    ap_assert(!is.null(meta), "This ordination carries no metadata, so `group` cannot be used.")
    ap_assert(group %in% names(meta), "Variable `{group}` is not in the sample metadata.")
    df$grp <- factor(as.character(meta[[group]][match(df$sample_id, rownames(meta))]))
    df <- df[!is.na(df$grp), , drop = FALSE]
  }
  if (!is.null(shape)) {
    ap_assert(shape %in% names(meta), "Variable `{shape}` is not in the sample metadata.")
    df$shp <- factor(as.character(meta[[shape]][match(df$sample_id, rownames(meta))]))
  }

  if (publication) {
    return(ap_plot_ordination_pub(df, ord, group, permanova, ellipse, axes, pub))
  }

  aes_args <- list(x = quote(.data$x), y = quote(.data$y))
  if (!is.null(group)) aes_args$colour <- quote(.data$grp)
  if (!is.null(shape)) aes_args$shape <- quote(.data$shp)

  p <- ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_point(size = point_size, alpha = 0.75)

  if (ellipse && !is.null(group)) {
    p <- p + ggplot2::stat_ellipse(level = 0.95, linewidth = 0.5, type = "t")
  }
  if (!is.null(group)) p <- p + ap_scale_colour(nlevels(df$grp), name = group)

  p <- p +
    ggplot2::labs(
      x = ap_axis_label(ord, axes[1]),
      y = ap_axis_label(ord, axes[2]),
      subtitle = ap_ord_subtitle(ord),
      caption = ap_ord_caption(ord)
    ) +
    ggplot2::coord_fixed() +
    ap_theme(grid = "both")

  if (!is.null(permanova) && !is.null(group)) {
    label <- ap_permanova_label(permanova, ord$metric, group)
    if (!is.null(label)) {
      # Headroom first, then the label. Without the expansion the annotation
      # lands on top of the ellipses whenever the cloud reaches the panel top,
      # and under coord_fixed it usually does.
      p <- p +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.28))) +
        ggplot2::annotate("label", x = -Inf, y = Inf, label = label,
                          hjust = -0.03, vjust = 1.05, size = 2.4,
                          fill = "white", alpha = 0.85, lineheight = 0.95)
    }
  }
  p
}

# The publication ordination. The statistics stay on the panel because a reader
# looks for them there; the verdict moves to the legend text with its reasoning,
# since a bare "confounded by dispersion" on a panel says too little to act on.
#' @keywords internal
ap_plot_ordination_pub <- function(df, ord, group, permanova, ellipse, axes, pub) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  if (is.null(group)) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(size = 0.9, alpha = 0.8, stroke = 0)
  } else {
    df$grp <- ap_pub_levels(df$grp, pub)
    cols <- ap_pub_palette(nlevels(df$grp), pub$palette)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y,
                                          colour = .data$grp, fill = .data$grp))
    if (ellipse) {
      p <- p +
        ggplot2::stat_ellipse(geom = "polygon", level = 0.95, type = "t",
                              alpha = 0.08, colour = NA) +
        ggplot2::stat_ellipse(level = 0.95, type = "t", linewidth = 0.4)
    }
    p <- p +
      ggplot2::geom_point(size = 0.9, alpha = 0.8, stroke = 0) +
      ggplot2::scale_colour_manual(values = cols, name = ap_pub_label(group, pub)) +
      ggplot2::scale_fill_manual(values = cols, name = ap_pub_label(group, pub))
  }
  p <- p +
    ggplot2::labs(x = ap_axis_label(ord, axes[1]), y = ap_axis_label(ord, axes[2])) +
    ggplot2::coord_fixed() +
    ap_theme_pub(legend = "right")

  if (!is.null(permanova) && !is.null(group)) {
    lines <- ap_permanova_plotmath(permanova, ord$metric, group)
    if (!is.null(lines)) {
      p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.25)))
      # One grob per line, since plotmath has no line break. Each baseline sits a fixed
      # 2.8 mm below the last; stepping by vjust instead spaces them by each line's own
      # height, which the superscript makes uneven.
      for (i in seq_along(lines)) {
        grob <- grid::textGrob(
          parse(text = lines[i])[[1]],
          x = grid::unit(1.5, "mm"), y = grid::unit(1, "npc") - grid::unit(3 + 2.8 * (i - 1), "mm"),
          hjust = 0, vjust = 0,
          gp = grid::gpar(fontsize = ap_pub_stat_pt, fontfamily = "Arial", col = "black"))
        p <- p + ggplot2::annotation_custom(grob)
      }
    }
  }
  attr(p, "ap_pub_size") <- ap_pub_ord_size(p)
  p
}

# Size for a coord_fixed ordination: one column (Bharat, 2026-09-18), with a
# height from the built panel's own x and y ranges (ellipses and headroom
# included), so the fixed aspect does not leave the page empty. Measured on the
# Baxter figure: the panel gets about 54 mm of the 89 mm width, and the x axis
# and margins take about 15 mm of height.
#' @keywords internal
ap_pub_ord_size <- function(p) {
  pp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  h <- 54 * diff(pp$y.range) / diff(pp$x.range) + 15
  list(width = "single", height = round(min(247, max(50, h))))
}

# In-panel statistics at 6 pt on two lines (Bharat's choice, 2026-09-18, after
# seeing 7 pt on three): at one column the panel is about 54 mm wide, and the
# PERMANOVA line needs about 64 mm at 7 pt. Nature allows 5 pt; Elsevier asks for
# 7 pt text, with 6 pt only for sub- and superscripts.
ap_pub_stat_pt <- 6

# PERMANOVA and dispersion as two plotmath expressions, so R-squared is set with a
# real superscript and the statistic symbols in italics.
#' @keywords internal
ap_permanova_plotmath <- function(pn, metric, term) {
  r <- pn$results[pn$results$metric == metric & pn$results$term == term, ]
  if (nrow(r) == 0L) return(NULL)
  d <- pn$dispersion[pn$dispersion$metric == metric & pn$dispersion$term == term, ]
  pv <- function(p) sprintf("italic(p) == '%s'", format.pval(p, digits = 2))
  c(sprintf("'PERMANOVA:' ~ italic(R)^2 == '%.4f' * ',' ~ italic(F) == '%.2f' * ',' ~ %s",
            r$R2[1], r$pseudo_F[1], pv(r$p[1])),
    if (nrow(d) == 0L || is.na(d$dispersion_p[1])) "'betadisper:' ~ '-'"
    else sprintf("'betadisper:' ~ %s", pv(d$dispersion_p[1])))
}

#' Legend text for a publication ordination
#'
#' What [ap_plot_ordination()] leaves off the publication panel: the method and
#' variance explained, the projection caveat when there is one, and the
#' PERMANOVA and dispersion results with the verdict and its explanation.
#'
#' @param ord An `ap_ordination` from [ap_ordinate()].
#' @param permanova An `ap_permanova` result, or `NULL`.
#' @param group The term tested.
#' @param pub [ap_pub_options()], for the group names in the pairwise lines.
#' @return A character vector, one line per entry.
#' @export
ap_ordination_legend <- function(ord, permanova = NULL, group = NULL, pub = ap_pub_options()) {
  out <- c(paste0(ap_ord_subtitle(ord), "."), ap_ord_caption(ord))
  if (!is.null(permanova) && !is.null(group)) {
    r <- permanova$results[permanova$results$metric == ord$metric &
                             permanova$results$term == group, ]
    v <- permanova$interpretation[permanova$interpretation$metric == ord$metric &
                                    permanova$interpretation$term == group, ]
    if (nrow(r) > 0L) {
      out <- c(out, sprintf("PERMANOVA (%s permutations): R2 = %.4f, pseudo-F = %.2f, p = %s, adjusted p = %s.",
                            format(permanova$permutations, big.mark = ","), r$R2[1],
                            r$pseudo_F[1], format.pval(r$p[1], digits = 2),
                            format.pval(r$p_adj[1], digits = 2)))
    }
    if (nrow(v) > 0L) {
      out <- c(out, paste0("Verdict: ", v$verdict[1], ". ", v$interpretation[1]))
    }
    out <- c(out, ap_pairwise_lines(permanova, ord$metric, group, pub))
  }
  out
}

# One legend line per pairwise PERMANOVA comparison for a metric and term, or
# nothing when the term has fewer than three groups (no pairwise was run).
#' @keywords internal
ap_pairwise_lines <- function(permanova, metric, term, pub = NULL) {
  pw <- permanova$pairwise
  if (is.null(pw)) return(NULL)
  pw <- pw[pw$metric == metric & pw$term == term, , drop = FALSE]
  if (nrow(pw) == 0L) return(NULL)
  c("Pairwise PERMANOVA (FDR-adjusted within this comparison):",
    sprintf("  %s vs %s: R2 = %.4f, pseudo-F = %.2f, p = %s, adjusted p = %s.",
            ap_pub_level_text(pw$group1, pub), ap_pub_level_text(pw$group2, pub), pw$R2, pw$pseudo_F,
            format.pval(pw$p, digits = 2), format.pval(pw$p_adj, digits = 2)))
}

#' @keywords internal
ap_ord_subtitle <- function(ord) {
  if (ord$method == "pcoa") {
    sprintf("PCoA on %s, first two axes explain %.1f%%",
            ord$metric, 100 * sum(ord$prop_explained[1:2]))
  } else {
    sprintf("NMDS on %s, stress = %.3f", ord$metric, ord$stress)
  }
}

#' @keywords internal
ap_ord_caption <- function(ord) {
  if (ord$method == "pcoa" && !is.na(ord$negative_eigenvalue_fraction) &&
      ord$negative_eigenvalue_fraction > 0.05) {
    sprintf("%.1f%% negative eigenvalue mass: non-Euclidean, projection distorts.",
            100 * ord$negative_eigenvalue_fraction)
  } else if (ord$method == "nmds") {
    "Stress < 0.1 good, 0.1-0.2 usable, > 0.2 not a map."
  } else {
    NULL
  }
}

#' @keywords internal
ap_permanova_label <- function(pn, metric, term) {
  r <- pn$results[pn$results$metric == metric & pn$results$term == term, ]
  if (nrow(r) == 0L) return(NULL)
  v <- pn$interpretation[pn$interpretation$metric == metric &
                           pn$interpretation$term == term, ]
  d <- pn$dispersion[pn$dispersion$metric == metric & pn$dispersion$term == term, ]
  sprintf("PERMANOVA: R2 = %.4f, F = %.2f, p = %s\nbetadisper p = %s\n%s",
          r$R2[1], r$pseudo_F[1], format.pval(r$p[1], digits = 2),
          if (nrow(d) == 0L || is.na(d$dispersion_p[1])) "-" else
            format.pval(d$dispersion_p[1], digits = 2),
          if (nrow(v) == 0L) "" else toupper(v$verdict[1]))
}

#' Plot within-group dispersion
#'
#' Distance from each sample to its group centroid, which is what `betadisper`
#' tests. This is the figure that shows whether a significant PERMANOVA is a
#' shift in location or a difference in spread.
#'
#' @param permanova An `ap_permanova` result from [ap_permanova()].
#' @param metric Which metric to plot. Defaults to the first.
#' @param term Which term to plot. Defaults to the first.
#' @param publication Draw the publication version: no subtitle or caption,
#'   [ap_theme_pub()], ggsci colours, and the betadisper statistics on the panel
#'   as plotmath. The method note and the PERMANOVA verdict go to the legend text
#'   from [ap_dispersion_legend()]. Default `FALSE`.
#' @param pub [ap_pub_options()] for `publication = TRUE`.
#'
#' @return A ggplot object.
#' @export
ap_plot_dispersion <- function(permanova, metric = NULL, term = NULL,
                               publication = FALSE, pub = ap_pub_options()) {
  ap_assert(inherits(permanova, "ap_permanova"),
            "`permanova` must come from `ap_permanova()`, not {class(permanova)[1]}.")
  metric <- metric %||% permanova$metrics[1]
  term <- term %||% permanova$terms[1]

  beta <- permanova$beta
  ap_assert(metric %in% names(beta$distances),
            "Metric `{metric}` is not in the underlying ap_beta object.")
  d <- beta$distances[[metric]]
  ids <- attr(d, "Labels")
  g <- beta$metadata[[term]][match(ids, rownames(beta$metadata))]
  ap_assert(!is.numeric(g),
            "`{term}` is continuous. Dispersion is defined between groups, not along a gradient.")

  keep <- !is.na(g)
  d <- stats::as.dist(as.matrix(d)[keep, keep])
  g <- droplevels(factor(as.character(g[keep])))

  set.seed(permanova$seed)
  bd <- vegan::betadisper(d, g)
  df <- data.frame(distance = bd$distances, grp = bd$group, stringsAsFactors = FALSE)

  if (publication) {
    return(ap_plot_dispersion_pub(df, permanova, metric, term, pub))
  }

  disp <- permanova$dispersion[permanova$dispersion$metric == metric &
                                 permanova$dispersion$term == term, ]
  label <- if (nrow(disp) > 0L && !is.na(disp$dispersion_p[1])) {
    sprintf("betadisper F = %.3f, p = %s\nspread ratio %.2fx",
            disp$dispersion_F[1], format.pval(disp$dispersion_p[1], digits = 2),
            disp$max_centroid_ratio[1])
  } else NULL

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$grp, y = .data$distance,
                                        fill = .data$grp)) +
    ggplot2::geom_boxplot(alpha = 0.55, linewidth = 0.3, colour = "grey25",
                          outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.16, height = 0, size = 0.5, alpha = 0.35,
                         colour = "grey15") +
    ap_scale_fill(nlevels(df$grp), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.22))) +
    ggplot2::labs(x = term, y = "Distance to group centroid",
                  subtitle = paste0("Within-group dispersion, ", metric),
                  caption = "PERMANOVA assumes these are comparable.") +
    ap_theme()

  if (!is.null(label)) {
    p <- p + ggplot2::annotate("text", x = -Inf, y = Inf, label = label,
                               hjust = -0.08, vjust = 1.2, size = 2.4,
                               colour = "grey20", lineheight = 0.95)
  }
  p
}

# The publication dispersion. Same treatment as the publication ordination: the
# statistics stay on the panel where the reader looks for them, and the method
# note and PERMANOVA verdict move to the legend text (ap_dispersion_legend()).
#' @keywords internal
ap_plot_dispersion_pub <- function(df, permanova, metric, term, pub) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  df$grp <- ap_pub_levels(df$grp, pub)
  cols <- ap_pub_palette(nlevels(df$grp), pub$palette)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$grp, y = .data$distance,
                                        fill = .data$grp)) +
    ggplot2::geom_boxplot(alpha = 0.6, linewidth = 0.3, colour = "grey25",
                          outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 0.6, alpha = 0.5,
                         colour = "black", stroke = 0) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30))) +
    ggplot2::labs(x = ap_pub_label(term, pub), y = "Distance to group centroid") +
    ap_theme_pub()

  lines <- ap_dispersion_plotmath(permanova, metric, term)
  if (!is.null(lines)) {
    # One grob per line, since plotmath has no line break. Baselines a fixed
    # 2.8 mm apart, as on the publication ordination.
    for (i in seq_along(lines)) {
      grob <- grid::textGrob(
        parse(text = lines[i])[[1]],
        x = grid::unit(1.5, "mm"), y = grid::unit(1, "npc") - grid::unit(3 + 2.8 * (i - 1), "mm"),
        hjust = 0, vjust = 0,
        gp = grid::gpar(fontsize = ap_pub_stat_pt, fontfamily = "Arial", col = "black"))
      p <- p + ggplot2::annotation_custom(grob)
    }
  }
  # A single boxplot panel at one column, the height fixed like the alpha panels.
  attr(p, "ap_pub_size") <- list(width = "single", height = 60)
  p
}

# The two tests the dispersion figure exists to compare, plus the spread ratio,
# as three plotmath lines: PERMANOVA on top (a location shift), betadisper below
# (a spread difference), so a reader sees both before reading the boxes. The
# statistic symbols are italic and the multiplication sign is a real glyph.
#' @keywords internal
ap_dispersion_plotmath <- function(pn, metric, term) {
  d <- pn$dispersion[pn$dispersion$metric == metric & pn$dispersion$term == term, ]
  if (nrow(d) == 0L || is.na(d$dispersion_p[1])) return(NULL)
  pv <- function(p) format.pval(p, digits = 2)
  r <- pn$results[pn$results$metric == metric & pn$results$term == term, ]
  perm <- if (nrow(r) > 0L) {
    sprintf("'PERMANOVA:' ~ italic(R)^2 == '%.4f' * ',' ~ italic(F) == '%.2f' * ',' ~ italic(p) == '%s'",
            r$R2[1], r$pseudo_F[1], pv(r$p[1]))
  }
  c(perm,
    sprintf("'betadisper:' ~ italic(F) == '%.2f' * ',' ~ italic(p) == '%s'",
            d$dispersion_F[1], pv(d$dispersion_p[1])),
    sprintf("'spread ratio' ~ '%.2f×'", d$max_centroid_ratio[1]))
}

#' Legend text for a publication dispersion figure
#'
#' What [ap_plot_dispersion()] leaves off the publication panel: what the
#' distances are, the note that PERMANOVA assumes equal dispersion, and the
#' betadisper result with the PERMANOVA verdict and its reasoning.
#'
#' @param permanova An `ap_permanova` result from [ap_permanova()].
#' @param metric Which metric. Defaults to the first.
#' @param term Which term. Defaults to the first.
#' @param pub [ap_pub_options()], for the group names in the pairwise lines.
#' @return A character vector, one line per entry.
#' @export
ap_dispersion_legend <- function(permanova, metric = NULL, term = NULL, pub = ap_pub_options()) {
  ap_assert(inherits(permanova, "ap_permanova"),
            "`permanova` must come from `ap_permanova()`, not {class(permanova)[1]}.")
  metric <- metric %||% permanova$metrics[1]
  term <- term %||% permanova$terms[1]
  out <- c(sprintf("Distance from each sample to its group centroid on %s, the quantity betadisper tests.",
                   metric),
           "PERMANOVA assumes groups have equal dispersion, so a location shift and a spread difference are read together unless this figure separates them.")
  d <- permanova$dispersion[permanova$dispersion$metric == metric &
                              permanova$dispersion$term == term, ]
  if (nrow(d) > 0L && !is.na(d$dispersion_p[1])) {
    out <- c(out, sprintf("betadisper: F = %.3f, p = %s, spread differs by %.2fx between groups.",
                          d$dispersion_F[1], format.pval(d$dispersion_p[1], digits = 2),
                          d$max_centroid_ratio[1]))
  }
  v <- permanova$interpretation[permanova$interpretation$metric == metric &
                                  permanova$interpretation$term == term, ]
  if (nrow(v) > 0L) {
    out <- c(out, paste0("Verdict: ", v$verdict[1], ". ", v$interpretation[1]))
  }
  c(out, ap_pairwise_lines(permanova, metric, term, pub))
}
