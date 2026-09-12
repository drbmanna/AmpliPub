# Figure defaults.
#
# One theme across every AmpliPub plot, so a figure panel assembled from
# several of them reads as one figure rather than three. Sized for a journal
# column at the default 3.5 inch width: a base size of 9 pt survives the
# reduction most publishers apply, where the ggplot2 default of 11 does not.

#' AmpliPub plot theme
#'
#' A clean theme for figures headed to a manuscript. Applied by every `ap_plot_*`
#' function; use it directly on your own ggplots to match.
#'
#' @param base_size Base font size in points. Default `9`.
#' @param base_family Font family. Default `""`, the device default.
#' @param grid Which grid lines to keep: `"y"`, `"x"`, `"both"` or `"none"`.
#' @return A ggplot2 theme.
#' @export
ap_theme <- function(base_size = 9, base_family = "", grid = c("y", "x", "both", "none")) {
  grid <- match.arg(grid)
  th <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "grey20", linewidth = 0.4),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey20",
                                               linewidth = 0.4),
      strip.text = ggplot2::element_text(size = base_size * 0.95, margin = ggplot2::margin(3, 3, 3, 3)),
      axis.text = ggplot2::element_text(colour = "grey20"),
      legend.key.size = ggplot2::unit(0.8, "lines"),
      legend.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size * 1.1, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.9, colour = "grey30"),
      plot.caption = ggplot2::element_text(size = base_size * 0.8, colour = "grey30",
                                           hjust = 0)
    )
  switch(
    grid,
    y = th + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()),
    x = th + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()),
    both = th,
    none = th + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  )
}

# Okabe-Ito. Chosen because it stays distinguishable under the three common
# forms of colour vision deficiency and survives greyscale printing, which a
# rainbow palette does not.
#' @keywords internal
ap_palette <- function(n) {
  okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                 "#E69F00", "#56B4E9", "#F0E442", "#000000")
  if (n <= length(okabe_ito)) return(okabe_ito[seq_len(n)])
  grDevices::colorRampPalette(okabe_ito)(n)
}

#' @keywords internal
ap_scale_fill <- function(n, ...) {
  ggplot2::scale_fill_manual(values = ap_palette(n), ...)
}

#' @keywords internal
ap_scale_colour <- function(n, ...) {
  ggplot2::scale_colour_manual(values = ap_palette(n), ...)
}

# Metric labels that say what the number means. "q0" on an axis is a label only
# the author understands.
#' @keywords internal
ap_metric_label <- function(metric) {
  labels <- c(
    q0 = "Hill q0\n(observed richness)",
    q1 = "Hill q1\n(effective no. of species)",
    q2 = "Hill q2\n(dominant-weighted)",
    evenness = "Pielou's evenness",
    faith_pd = "Faith's PD",
    shannon_entropy = "Shannon entropy (bits)"
  )
  out <- unname(labels[metric])
  ifelse(is.na(out), metric, out)
}
