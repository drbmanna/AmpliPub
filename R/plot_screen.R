#' Plot an exploratory screen
#'
#' Effect size for the top rows of an [ap_screen()] result, with the stability
#' of each rank written next to it. The x axis is variance explained adjusted
#' for degrees of freedom, so the figure carries the ranking the screen uses and
#' never a p-value ranking.
#'
#' @param screen An `ap_screen` object from [ap_screen()].
#' @param n Number of top rows to show. Default `20`.
#' @param publication Draw the publication version instead: a heatmap for one `family`,
#'   one row per variable and one column per metric, so every test is shown and none is
#'   cut off. Fill is the adjusted effect size; a dot marks q < 0.05, an open circle a beta
#'   test that is significant but whose groups also differ in dispersion. Alpha and beta
#'   are separate figures, with their own colour scales, because the two effect sizes are
#'   proportions of different variances. Both list the variables in the same order.
#'   Counts, rank stability and caveats go to [ap_screen_legend()].
#' @param family For `publication = TRUE`: `"alpha"` or `"beta"`.
#' @param pub [ap_pub_options()] for `publication = TRUE`.
#'
#' @return A ggplot object.
#' @export
ap_plot_screen <- function(screen, n = 20L, publication = FALSE,
                           family = c("alpha", "beta"), pub = ap_pub_options()) {
  ap_assert(inherits(screen, "ap_screen"),
            "`screen` must come from `ap_screen()`, not {class(screen)[1]}.")
  if (publication) return(ap_plot_screen_pub(screen, match.arg(family), pub))
  r <- utils::head(screen$results, n)
  r$label <- paste0(r$variable, " | ", r$metric)
  r$label <- factor(r$label, levels = rev(unique(r$label)))
  r$family_label <- ifelse(r$family == "alpha",
                           "Alpha (eta2[H] or rho2)", "Beta (PERMANOVA R2)")
  r$stability_label <- sprintf("%.0f%%", 100 * r$stability)

  ggplot2::ggplot(r, ggplot2::aes(x = .data$effect_adj, y = .data$label,
                                  colour = .data$family_label)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$effect_adj, yend = .data$label),
                          x = 0, linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(alpha = .data$stability), size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$stability_label), hjust = -0.35,
                       size = 2.3, colour = "grey30", show.legend = FALSE) +
    ggplot2::scale_alpha_continuous(range = c(0.3, 1), limits = c(0, 1),
                                    name = sprintf("Top-%d stability", screen$top_k)) +
    ap_scale_colour(2, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.2))) +
    ggplot2::labs(
      x = "Variance explained, adjusted for degrees of freedom",
      y = NULL,
      subtitle = sprintf("Top %d of %d tests, ranked by adjusted effect size", nrow(r), screen$n_tests),
      caption = sprintf(paste0("Labels: share of %d subsamples (%.0f%% of samples) in the top %d. ",
                               "BH across all %d tests. Hypothesis-generating."),
                        screen$n_resample, 100 * screen$fraction, screen$top_k, screen$n_tests)
    ) +
    ap_theme(grid = "x")
}

# Publication heatmaps ----------------------------------------------------------------------------

#' @keywords internal
ap_screen_metric_labels <- function() {
  m <- c("q0", "q1", "q2", "evenness", "faith_pd", "bray_curtis", "jaccard",
         "unweighted_unifrac", "weighted_unifrac")
  stats::setNames(ap_metric_label_pub(m), m)
}

# Every variable the screen tested, in the screen's own order (strongest single test first).
# Both family figures use this order, so their rows line up side by side.
#' @keywords internal
ap_screen_variable_order <- function(screen) {
  r <- screen$results
  vars <- unique(screen$by_variable$variable %||% r$variable[order(-r$effect_adj)])
  c(vars, setdiff(unique(r$variable), vars))
}

#' @keywords internal
ap_plot_screen_pub <- function(screen, family, pub) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  vars <- ap_screen_variable_order(screen)
  r <- screen$results[screen$results$family == family, , drop = FALSE]
  ap_assert(nrow(r) > 0L, "The screen ran no {family} tests, so there is nothing to draw.")
  # Columns in a fixed order, so figures from different studies line up.
  known <- names(ap_screen_metric_labels())
  mets <- c(intersect(known, unique(r$metric)), setdiff(unique(r$metric), known))
  mlab <- ap_screen_metric_labels()[mets]
  mlab[is.na(mlab)] <- mets[is.na(mlab)]

  vlab <- vapply(vars, ap_pub_label, character(1), pub = pub)
  ap_assert(!anyDuplicated(vlab), "Two variables share a label in `labels`; rows would merge.")
  r$row <- factor(vlab[match(r$variable, vars)], levels = rev(vlab))
  r$col <- factor(unname(mlab[match(r$metric, mets)]), levels = unname(mlab))
  # Adjusted effect sizes can fall below zero; there is no negative variance to show.
  r$fill <- pmax(r$effect_adj, 0)
  sig <- !is.na(r$q) & r$q < 0.05
  conf <- sig & !is.na(r$verdict) & r$verdict == "confounded by dispersion"
  r$mark <- ifelse(conf, "Dispersion also differs", ifelse(sig, "q < 0.05", NA))
  marks <- r[!is.na(r$mark), , drop = FALSE]

  p <- ggplot2::ggplot(r, ggplot2::aes(x = .data$col, y = .data$row)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white", linewidth = 0.3) +
    ggplot2::geom_point(data = marks, ggplot2::aes(shape = .data$mark), size = 0.9,
                        stroke = 0.35, colour = "black") +
    ggplot2::scale_shape_manual(values = c(`q < 0.05` = 16, `Dispersion also differs` = 1),
                                breaks = intersect(c("q < 0.05", "Dispersion also differs"),
                                                   marks$mark),
                                name = NULL, drop = TRUE) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#2171B5", limits = c(0, NA),
                                 name = "Variance explained\n(adjusted)") +
    ggplot2::scale_x_discrete(expand = c(0, 0), position = "top") +
    ggplot2::scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL) +
    ap_theme_pub(legend = "bottom") +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x.top = ggplot2::element_text(size = 6, angle = 45, hjust = 0,
                                                           vjust = 0),
                   axis.text.y = ggplot2::element_text(size = 6),
                   legend.box = "vertical", legend.spacing.y = ggplot2::unit(1, "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.text = ggplot2::element_text(size = 6),
                   legend.title = ggplot2::element_text(size = 6),
                   legend.key.height = ggplot2::unit(2.5, "mm"),
                   legend.key.width = ggplot2::unit(6, "mm"),
                   # Room on the right for the last rotated column label.
                   plot.margin = ggplot2::margin(2, 12, 2, 2, unit = "mm"))
  attr(p, "ap_pub_size") <- list(width = "single",
                                 height = min(247, round(3.4 * length(vars) + 45)))
  p
}

#' Legend text for a publication screen figure
#'
#' What a screen heatmap leaves off: how many tests, the correction, what the effect size
#' is, the tests whose rank was stable under resampling, for beta diversity the tests that
#' differ in dispersion rather than (or as well as) location, and the variables that could
#' not be tested.
#'
#' @param screen An `ap_screen` object from [ap_screen()].
#' @param family `"alpha"` or `"beta"`, as drawn by `ap_plot_screen(publication = TRUE)`.
#' @param pub [ap_pub_options()], for variable labels.
#' @param min_stability Rank stability at or above which a test is named. Default `0.5`.
#' @return A character vector, one line per entry.
#' @export
ap_screen_legend <- function(screen, family = c("alpha", "beta"), pub = ap_pub_options(),
                             min_stability = 0.5) {
  ap_assert(inherits(screen, "ap_screen"),
            "`screen` must come from `ap_screen()`, not {class(screen)[1]}.")
  family <- match.arg(family)
  all <- screen$results
  r <- all[all$family == family, , drop = FALSE]
  lab <- function(v, m) sprintf("%s (%s)", vapply(v, ap_pub_label, character(1), pub = pub),
                                ap_screen_metric_labels()[m])
  stable <- r[!is.na(r$stability) & r$stability >= min_stability, , drop = FALSE]
  stable <- stable[order(-stable$stability), , drop = FALSE]
  conf <- r[!is.na(r$verdict) & r$verdict == "confounded by dispersion" &
              !is.na(r$q) & r$q < 0.05, , drop = FALSE]
  donly <- r[!is.na(r$verdict) & r$verdict == "dispersion only", , drop = FALSE]
  skipped <- screen$skipped
  other <- if (family == "alpha") "beta" else "alpha"
  c(
    sprintf(paste0("Exploratory screen, %s diversity: %d tests (%d variables by %d metrics). ",
                   "Benjamini-Hochberg correction across all %d tests of the screen, alpha and ",
                   "beta together; %d of these with q < 0.05 (dots)."),
            family, nrow(r), length(unique(r$variable)), length(unique(r$metric)), nrow(all),
            sum(!is.na(r$q) & r$q < 0.05)),
    if (family == "alpha") {
      c(paste0("Fill: variance explained, adjusted for degrees of freedom: eta squared from the ",
               "Kruskal-Wallis H for categorical variables, adjusted Spearman rho squared for ",
               "numeric ones. Values below zero are shown as zero."),
        ap_alpha_names_note(unique(r$metric)))
    } else {
      paste0("Fill: PERMANOVA R2 adjusted for degrees of freedom. Values below zero are shown ",
             "as zero.")
    },
    sprintf(paste0("The colour scale is not shared with the %s diversity figure: the two ",
                   "measure proportions of different variances."), other),
    if (nrow(stable)) {
      sprintf(paste0("Stable in rank (in the top %d %s tests in at least %.0f%% of %d ",
                     "subsamples of %.0f%% of samples): %s."),
              as.integer(screen$top_k), family, 100 * min_stability,
              as.integer(screen$n_resample), 100 * screen$fraction,
              paste(sprintf("%s %.0f%%", lab(stable$variable, stable$metric),
                            100 * stable$stability), collapse = "; "))
    } else if (all(is.na(r$stability))) {
      sprintf(paste0("Rank stability not estimated: with %d %s tests, every test is in the ",
                     "top %d by construction."), nrow(r), family, as.integer(screen$top_k))
    } else {
      sprintf(paste0("No %s test stayed in the top %d %s tests in at least %.0f%% of ",
                     "subsamples."), family, as.integer(screen$top_k), family,
              100 * min_stability)
    },
    if (nrow(conf)) {
      paste0("Open circles: significant, but group dispersions also differ (betadisper ",
             "p < 0.05), so location and spread cannot be separated: ",
             paste(lab(conf$variable, conf$metric), collapse = "; "), ".")
    },
    if (nrow(donly)) {
      paste0("Dispersion differs without a location shift: ",
             paste(lab(donly$variable, donly$metric), collapse = "; "), ".")
    },
    if (!is.null(skipped) && nrow(skipped)) {
      paste0("Not tested: ", paste(sprintf("%s (%s)", skipped$variable, skipped$reason),
                                   collapse = "; "), ".")
    },
    "Hypothesis-generating: these tests were not planned, and are not reported as findings."
  )
}
