#' Plot an exploratory screen
#'
#' Effect size for the top rows of an [ap_screen()] result, with the stability
#' of each rank written next to it. The x axis is variance explained adjusted
#' for degrees of freedom, so the figure carries the ranking the screen uses and
#' never a p-value ranking.
#'
#' @param screen An `ap_screen` object from [ap_screen()].
#' @param n Number of top rows to show. Default `20`.
#'
#' @return A ggplot object.
#' @export
ap_plot_screen <- function(screen, n = 20L) {
  ap_assert(inherits(screen, "ap_screen"),
            "`screen` must come from `ap_screen()`, not {class(screen)[1]}.")
  r <- utils::head(screen$results, n)
  r$label <- paste0(r$variable, " | ", r$metric)
  r$label <- factor(r$label, levels = rev(unique(r$label)))
  r$family_label <- ifelse(r$family == "alpha",
                           "Alpha (epsilon2 or rho2)", "Beta (PERMANOVA R2)")
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
