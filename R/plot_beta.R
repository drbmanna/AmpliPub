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
#'
#' @return A ggplot object.
#' @export
ap_plot_ordination <- function(ord, group = NULL, shape = NULL, permanova = NULL,
                               ellipse = TRUE, axes = c(1, 2), point_size = 1.6) {
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
#'
#' @return A ggplot object.
#' @export
ap_plot_dispersion <- function(permanova, metric = NULL, term = NULL) {
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
