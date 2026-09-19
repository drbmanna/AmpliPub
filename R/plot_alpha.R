#' Plot alpha diversity by group
#'
#' One panel per metric, with the test result and effect size written onto each
#' panel rather than left to the caption. A reviewer should be able to read the
#' claim and its strength off the figure.
#'
#' The effect size is shown alongside the p-value on purpose. A boxplot with
#' `p < 0.001` and nothing else invites the reader to assume a large difference,
#' which at n in the hundreds it usually is not.
#'
#' @param alpha An `ap_alpha` object from [ap_alpha()].
#' @param group Metadata variable to split by.
#' @param test Optional `ap_alpha_test` result to annotate with. When `NULL`
#'   (default) and `annotate = TRUE`, [ap_alpha_test()] is run for you.
#' @param type `"box"` (default) or `"violin"`.
#' @param points Draw the individual samples. Default `TRUE`.
#' @param metrics Metrics to include. Defaults to everything in `alpha`.
#' @param annotate Write the test result onto each panel. Default `TRUE`.
#' @param ncol Facet columns. `NULL` lets ggplot2 decide.
#' @param publication Draw the publication version: no caption and no
#'   statistics on the panels, [ap_theme_pub()], a ggsci palette, each metric's
#'   name as its panel's y-axis title, and a panel grid chosen from the number
#'   of metrics. The statistics and the caveat go to the legend text from
#'   [ap_alpha_legend()]. Default `FALSE`.
#' @param pub [ap_pub_options()] for `publication = TRUE`: palette, axis
#'   titles for metadata variables, and level capitalization.
#'
#' @return A ggplot object.
#' @export
ap_plot_alpha <- function(alpha,
                          group,
                          test = NULL,
                          type = c("box", "violin"),
                          points = TRUE,
                          metrics = NULL,
                          annotate = TRUE,
                          ncol = NULL,
                          publication = FALSE,
                          pub = ap_pub_options()) {
  ap_assert(inherits(alpha, "ap_alpha"),
            "`alpha` must come from `ap_alpha()`, not {class(alpha)[1]}.")
  type <- match.arg(type)
  meta <- alpha$metadata
  ap_assert(group %in% names(meta),
            "Variable `{group}` is not in the sample metadata.")

  metrics <- metrics %||% unique(alpha$values$metric)
  df <- alpha$values[alpha$values$metric %in% metrics, ]
  df$grp <- meta[[group]][match(df$sample_id, rownames(meta))]
  df <- df[!is.na(df$grp) & !is.na(df$value), , drop = FALSE]
  ap_assert(nrow(df) > 0L, "Nothing to plot: every value is missing for `{group}`.")

  df$grp <- factor(as.character(df$grp))
  df$metric <- factor(df$metric, levels = metrics,
                      labels = if (publication) ap_metric_label_pub(metrics) else ap_metric_label(metrics))

  if (publication) return(ap_plot_alpha_pub(df, type, points, pub, group))

  if (annotate && is.null(test)) {
    test <- tryCatch(ap_alpha_test(alpha, group, metrics = metrics),
                     error = function(e) NULL)
  }

  n_grp <- nlevels(df$grp)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$grp, y = .data$value, fill = .data$grp))

  if (type == "violin") {
    # trim = TRUE. The kernel tails of an untrimmed violin run past the observed
    # range, which on richness or evenness draws density at values the quantity
    # cannot take, such as a negative number of species.
    p <- p + ggplot2::geom_violin(alpha = 0.55, linewidth = 0.3, colour = "grey25",
                                  scale = "width", trim = TRUE)
    p <- p + ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
                                   linewidth = 0.3)
  } else {
    p <- p + ggplot2::geom_boxplot(alpha = 0.55, linewidth = 0.3, colour = "grey25",
                                   outlier.shape = if (points) NA else 19,
                                   outlier.size = 0.6)
  }
  if (points) {
    p <- p + ggplot2::geom_jitter(width = 0.16, height = 0, size = 0.5,
                                  alpha = 0.35, colour = "grey15")
  }

  p <- p +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = ncol) +
    ap_scale_fill(n_grp, guide = "none") +
    ggplot2::labs(x = group, y = NULL,
                  caption = ap_alpha_caption(alpha)) +
    ap_theme()

  if (!is.null(test)) {
    # Annotation is pinned to the panel corner, and the top of the panel is
    # expanded to make room for it. Placing it at a data coordinate instead lets
    # it land on top of the distribution whenever the data reach the top of the
    # panel, which under free y scales is most of them.
    p <- p +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.30))) +
      ggplot2::geom_text(
        data = ap_alpha_annotations(test, df, metrics),
        mapping = ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
        inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 2.3,
        colour = "grey20", lineheight = 0.95
      )
  }
  p
}

# The publication version. Boxes are filled opaque: with a translucent fill the
# whisker segment ggplot2 draws behind the box shows through as a line across it.
#' @keywords internal
ap_plot_alpha_pub <- function(df, type, points, pub, group) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  df$grp <- ap_pub_levels(df$grp, pub)
  n_grp <- nlevels(df$grp)
  facet <- ap_pub_facet("metric", nlevels(droplevels(df$metric)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$grp, y = .data$value, fill = .data$grp))
  if (type == "violin") {
    p <- p +
      ggplot2::geom_violin(linewidth = 0.3, colour = "black", scale = "width", trim = TRUE) +
      ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = 0.3)
  } else {
    p <- p +
      ggplot2::stat_boxplot(geom = "errorbar", width = 0.25, linewidth = 0.3) +
      ggplot2::geom_boxplot(linewidth = 0.3, colour = "black", width = 0.6,
                            outlier.shape = if (points) NA else 19, outlier.size = 0.5)
  }
  if (points) {
    p <- p + ggplot2::geom_jitter(width = 0.15, height = 0, size = 0.6, alpha = 0.5,
                                  colour = "black", stroke = 0)
  }
  p <- p +
    facet[[1]] +
    ggplot2::scale_fill_manual(values = ap_pub_palette(n_grp, pub$palette), guide = "none") +
    ggplot2::labs(x = ap_pub_label(group, pub), y = NULL) +
    ap_theme_pub()
  attr(p, "ap_pub_size") <- facet[[2]][c("width", "height")]
  p
}

#' Legend text for a publication alpha-diversity figure
#'
#' Everything [ap_plot_alpha()] writes on the report figure and leaves off the
#' publication figure: the rarefaction statement or the caveat that the data
#' were not rarefied, and for each metric the test, n, adjusted p and effect
#' size with its interval. Written next to the figure by [ap_save_figure()].
#'
#' @param alpha An `ap_alpha` object.
#' @param test An `ap_alpha_test` result, or `NULL`.
#' @return A character vector, one line per entry.
#' @export
ap_alpha_legend <- function(alpha, test = NULL) {
  out <- c(ap_alpha_caption(alpha), ap_alpha_names_note(unique(alpha$values$metric)))
  if (!is.null(test)) {
    r <- test$results
    out <- c(out, vapply(seq_len(nrow(r)), function(i) {
      sprintf("%s: %s, n = %d, q = %s, %s = %.3f%s.",
              ap_metric_label_pub(r$metric[i]), r$test[i], r$n[i],
              format.pval(r$p_adj[i], digits = 2), r$effect[i], r$estimate[i],
              if (is.na(r$ci_low[i])) "" else
                sprintf(" (95%% CI %.3f to %.3f)", r$ci_low[i], r$ci_high[i]))
    }, character(1)))
  }
  out
}

# The publication figures say "Shannon" and "Inverse Simpson"; the legend says what number
# that is, and why Chao1 is absent when it is.
#' @keywords internal
ap_alpha_names_note <- function(metrics) {
  c(
    if (any(c("q1", "q2") %in% metrics)) {
      paste0("Shannon and inverse Simpson are the Hill numbers of order 1 and 2 (exp(H) and ",
             "1/sum(p^2)), in effective numbers of features; richness is the order-0 Hill number.")
    },
    if (!any(c("chao1", "ace") %in% metrics)) {
      paste0("Chao1 and ACE are not reported: both estimate unseen richness from singletons, ",
             "which denoising removes.")
    }
  )
}

#' @keywords internal
ap_alpha_caption <- function(alpha) {
  if (isTRUE(alpha$rarefied)) {
    paste0("Rarefied to ", format(alpha$depth, big.mark = ","), " reads, ",
           alpha$n_iter, " iterations averaged, seed ", alpha$seed, ".",
           if (length(alpha$dropped) > 0L)
             paste0(" ", length(alpha$dropped), " sample(s) below depth excluded.") else "")
  } else if (!is.na(alpha$common_depth %||% NA_real_)) {
    paste0("All samples at a common depth of ", format(alpha$common_depth, big.mark = ","),
           " reads; not rarefied again here.")
  } else {
    "Not rarefied; richness is not comparable across unequal library sizes."
  }
}

# Annotation sits at the top-left of each panel, in panel coordinates, so it
# lands correctly under free y scales.
#' @keywords internal
ap_alpha_annotations <- function(test, df, metrics) {
  r <- test$results
  rows <- lapply(metrics, function(m) {
    hit <- r[r$metric == m, ]
    if (nrow(hit) == 0L) return(NULL)
    sub <- df[df$metric == ap_metric_label(m), ]
    if (nrow(sub) == 0L) return(NULL)
    label <- sprintf("%s, q = %s\n%s = %.3f%s",
                     hit$test[1],
                     format.pval(hit$p_adj[1], digits = 2),
                     hit$effect[1], hit$estimate[1],
                     if (is.na(hit$ci_low[1])) "" else
                       sprintf("\n[%.3f, %.3f]", hit$ci_low[1], hit$ci_high[1]))
    data.frame(metric = ap_metric_label(m), label = label,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  out$metric <- factor(out$metric, levels = levels(df$metric))
  out
}

#' Plot a rarefaction curve
#'
#' Richness against sequencing depth, one line per sample or one ribbon per
#' group. This is the figure that justifies a rarefaction depth: where the
#' curves flatten, deeper sequencing stops finding new features.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param group Optional metadata variable to colour and summarise by.
#' @param depths Depths to evaluate. Defaults to 10 steps up to the median
#'   library size.
#' @param n_iter Iterations averaged at each depth. Default `5`.
#' @param seed Random seed, shown in the caption.
#' @param metric Metric to curve. Default `"q0"`.
#'
#' @return A ggplot object.
#' @export
ap_plot_rarefaction <- function(x, group = NULL, depths = NULL, n_iter = 5L,
                                seed = 1L, metric = "q0") {
  counts <- SummarizedExperiment::assay(x, "counts")
  lib <- colSums(counts)
  if (is.null(depths)) {
    depths <- unique(round(seq(max(1, min(lib) / 10), stats::median(lib), length.out = 10)))
  }

  set.seed(seed)
  rows <- lapply(depths, function(d) {
    keep <- names(lib)[lib >= d]
    if (length(keep) == 0L) return(NULL)
    sub_counts <- counts[, keep, drop = FALSE]
    acc <- NULL
    for (i in seq_len(n_iter)) {
      one <- ap_alpha_one(ap_rarefy_matrix(sub_counts, d), metric, NULL)
      acc <- if (is.null(acc)) one else {acc$value <- acc$value + one$value; acc}
    }
    acc$value <- acc$value / n_iter
    acc$depth <- d
    acc
  })
  df <- do.call(rbind, Filter(Negate(is.null), rows))

  meta <- as.data.frame(SummarizedExperiment::colData(x))
  if (!is.null(group)) {
    ap_assert(group %in% names(meta), "Variable `{group}` is not in the sample metadata.")
    df$grp <- factor(as.character(meta[[group]][match(df$sample_id, rownames(meta))]))
    df <- df[!is.na(df$grp), , drop = FALSE]
    summ <- stats::aggregate(value ~ depth + grp, data = df, FUN = function(v) {
      c(mean = mean(v), lo = stats::quantile(v, 0.25), hi = stats::quantile(v, 0.75))
    })
    summ <- data.frame(depth = summ$depth, grp = summ$grp,
                       as.data.frame(summ$value))
    names(summ)[3:5] <- c("mean", "lo", "hi")
    p <- ggplot2::ggplot(summ, ggplot2::aes(x = .data$depth, y = .data$mean,
                                            colour = .data$grp, fill = .data$grp)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                           alpha = 0.18, colour = NA) +
      ggplot2::geom_line(linewidth = 0.6) +
      ap_scale_colour(nlevels(summ$grp), name = group) +
      ap_scale_fill(nlevels(summ$grp), name = group)
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$depth, y = .data$value,
                                          group = .data$sample_id)) +
      ggplot2::geom_line(alpha = 0.2, linewidth = 0.25, colour = "#0072B2")
  }

  p +
    ggplot2::scale_x_continuous(labels = scales::comma) +
    ggplot2::labs(x = "Sequencing depth (reads)", y = ap_metric_label(metric),
                  caption = paste0("Mean of ", n_iter, " rarefactions per depth, seed ",
                                   seed, ". Samples below a depth are excluded at that depth, ",
                                   "so curves shorten rather than bend.")) +
    ap_theme(grid = "both")
}
