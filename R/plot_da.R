#' Volcano plot, faceted by method
#'
#' Effect size against adjusted significance, one panel per method.
#'
#' Panels do not share an x scale. The four methods report effects on four
#' different quantities (natural-log fold change, log2 fold change, a
#' coefficient on log relative abundance, a standardised CLR difference), and a
#' shared axis would present them as if one number meant the same thing in each.
#'
#' @param da An `ap_da` result from [ap_da()].
#' @param contrast Which contrast to plot. Defaults to the first.
#' @param label_top Number of features to label per panel. Default `8`.
#' @param label_by `"effect"` (default) labels the largest effects,
#'   `"significance"` the smallest q-values.
#'
#' @return A ggplot object.
#' @export
ap_plot_volcano <- function(da, contrast = NULL, label_top = 8L,
                            label_by = c("effect", "significance")) {
  ap_assert(inherits(da, "ap_da"), "`da` must come from `ap_da()`, not {class(da)[1]}.")
  label_by <- match.arg(label_by)

  df <- da$results
  contrast <- contrast %||% df$contrast[1]
  df <- df[df$contrast == contrast & !is.na(df$p_adj), , drop = FALSE]
  ap_assert(nrow(df) > 0L, "No results for contrast `{contrast}`.")

  df$neglog_q <- -log10(pmax(df$p_adj, .Machine$double.xmin))
  df$call <- ifelse(!df$significant, "not significant",
                    ifelse(df$effect > 0, "enriched", "depleted"))
  df$call <- factor(df$call, levels = c("depleted", "not significant", "enriched"))
  df$label <- df$taxon_label %||% df$feature

  labels <- do.call(rbind, lapply(split(df, df$method), function(d) {
    d <- d[d$significant, , drop = FALSE]
    if (nrow(d) == 0L) return(NULL)
    d <- d[order(if (label_by == "effect") -abs(d$effect) else d$p_adj), ]
    utils::head(d, label_top)
  }))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$effect, y = .data$neglog_q,
                                        colour = .data$call)) +
    ggplot2::geom_hline(yintercept = -log10(da$alpha), linetype = "dashed",
                        colour = "grey55", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
    ggplot2::geom_point(size = 0.9, alpha = 0.75) +
    ggplot2::facet_wrap(~ method, scales = "free_x") +
    ggplot2::scale_colour_manual(
      values = c(depleted = "#0072B2", `not significant` = "grey75", enriched = "#D55E00"),
      name = NULL, drop = FALSE
    ) +
    ggplot2::labs(
      x = "Effect size (each method on its own scale)",
      y = expression(-log[10] * " adjusted p"),
      subtitle = paste0(contrast, ", reference ", da$reference),
      caption = ap_wrap(paste0(
        "Panels do not share an x scale: the methods report different quantities. ",
        "Dashed line is the ", format(da$alpha), " threshold on the ",
        da$p_adj_method, "-adjusted p-value. Prevalence filter ", da$prv_cut,
        " applied once to all methods."))
    ) +
    ap_theme(grid = "both")

  if (!is.null(labels) && nrow(labels) > 0L &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = labels, mapping = ggplot2::aes(label = .data$label),
      size = 2, max.overlaps = 20, min.segment.length = 0.2,
      segment.size = 0.2, colour = "grey20", show.legend = FALSE
    )
  }
  p
}

#' Effect sizes with intervals for the consensus set
#'
#' Point estimate and 95% interval per method for the features the methods
#' agreed on, ordered by the size of the effect.
#'
#' This is the figure a differential abundance result should lead with. A list
#' of significant features with q-values says nothing about how large the
#' differences are, and at a few hundred samples most of them are small.
#'
#' @param concordance An `ap_da_concordance` result from [ap_da_concordance()].
#' @param features Features to show. Defaults to the consensus set.
#' @param max_features Cap on how many to draw. Default `25`.
#'
#' @return A ggplot object.
#' @export
ap_plot_da_effects <- function(concordance, features = NULL, max_features = 25L) {
  ap_assert(inherits(concordance, "ap_da_concordance"),
            "`concordance` must come from `ap_da_concordance()`, not {class(concordance)[1]}.")
  features <- features %||% concordance$consensus
  ap_assert(length(features) > 0L,
            paste0("The consensus set is empty, so there is nothing to plot. ",
                   "Lower `min_methods` in `ap_da_concordance()` to see partial agreement."))

  r <- concordance$da$results
  r <- r[r$contrast == concordance$contrast & r$feature %in% features, , drop = FALSE]

  order_by <- stats::aggregate(effect ~ feature, data = r, FUN = function(v) mean(abs(v)))
  order_by <- order_by[order(-order_by$effect), ]
  keep <- utils::head(order_by$feature, max_features)
  r <- r[r$feature %in% keep, , drop = FALSE]

  r$label <- r$taxon_label %||% r$feature
  label_order <- r$label[match(keep, r$feature)]
  r$label <- factor(r$label, levels = rev(label_order))
  r$lo <- r$effect - 1.96 * r$se
  r$hi <- r$effect + 1.96 * r$se

  ggplot2::ggplot(r, ggplot2::aes(x = .data$effect, y = .data$label,
                                  colour = .data$method)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    # geom_errorbarh() is deprecated in ggplot2 4.0; the horizontal form is now
    # geom_errorbar() with orientation = "y".
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                           orientation = "y", width = 0, linewidth = 0.4,
                           position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(size = 1.3, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::facet_wrap(~ method, scales = "free_x", nrow = 1) +
    ap_scale_colour(length(unique(r$method)), guide = "none") +
    ggplot2::labs(
      x = "Effect (each method on its own scale)", y = NULL,
      subtitle = paste0("Consensus features, ", concordance$contrast),
      caption = ap_wrap(paste0(
        "Bars are 95% intervals from each method's own standard error. ALDEx2's ",
        "is its within-condition dispersion, not a standard error, so its bar is ",
        "a spread rather than a confidence interval."))
    ) +
    ap_theme(grid = "x")
}

#' Concordance between methods
#'
#' Which methods called which features, as a presence grid ordered by how many
#' methods agreed.
#'
#' @param concordance An `ap_da_concordance` result from [ap_da_concordance()].
#' @param max_features Cap on features drawn. Default `40`.
#'
#' @return A ggplot object.
#' @export
ap_plot_concordance <- function(concordance, max_features = 40L) {
  ap_assert(inherits(concordance, "ap_da_concordance"),
            "`concordance` must come from `ap_da_concordance()`, not {class(concordance)[1]}.")
  feat <- concordance$features[concordance$features$n_methods > 0L, , drop = FALSE]
  ap_assert(nrow(feat) > 0L, "No feature was called by any method.")
  feat <- utils::head(feat, max_features)

  methods <- concordance$methods
  long <- do.call(rbind, lapply(methods, function(m) {
    data.frame(feature = feat$feature,
               label = feat$taxon_label %||% feat$feature,
               method = m,
               called = feat[[paste0("sig_", m)]],
               n_methods = feat$n_methods,
               stringsAsFactors = FALSE)
  }))
  long$label <- factor(long$label, levels = rev(unique(long$label[order(-long$n_methods)])))
  long$method <- factor(long$method, levels = methods)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$method, y = .data$label,
                                     fill = .data$called)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "grey92"),
                               name = NULL,
                               labels = c(`TRUE` = "called", `FALSE` = "not called")) +
    ggplot2::labs(
      x = NULL, y = NULL,
      subtitle = paste0("Called by method, ", concordance$contrast),
      caption = ap_wrap(paste0(
        length(concordance$consensus), " of ", nrow(concordance$features),
        " features were called by at least ", concordance$min_methods,
        " methods with a consistent direction. The rows called by one method ",
        "are the ones to be careful about."))
    ) +
    ap_theme(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
