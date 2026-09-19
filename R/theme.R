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
    shannon_entropy = "Shannon entropy (bits)",
    chao1 = "Chao1\n(estimated richness)",
    ace = "ACE\n(estimated richness)"
  )
  out <- unname(labels[metric])
  ifelse(is.na(out), metric, out)
}

# "Other" is a pooled residual, not a taxon, and it is often the largest band.
# Giving it a colour from the categorical palette makes it compete with real
# taxa for the reader's attention; giving it black, which is where Okabe-Ito
# ends, makes it dominate outright. It gets a light grey and is excluded from
# the palette so no real taxon is starved of a colour.
#' @keywords internal
ap_scale_fill_taxa <- function(levels, ...) {
  real <- setdiff(levels, "Other")
  values <- stats::setNames(ap_palette(length(real)), real)
  if ("Other" %in% levels) values <- c(values, Other = "#D9D9D9")
  ggplot2::scale_fill_manual(values = values, breaks = levels, ...)
}

# Publication figures ----------------------------------------------------------------
#
# The report figures above carry their statistics and caveats on the panel. The
# publication set does not: journals want the numbers in the figure legend, so
# every publication figure is saved with a `_legend.txt` holding what was taken
# off the panel (see ap_save_figure()). Sizes are for the figure drawn at its
# final printed width, so 7 pt here is 7 pt on paper.

ap_pub_require <- function(pkg) {
  ap_assert(requireNamespace(pkg, quietly = TRUE),
            "Publication figures need the {pkg} package. Install it, or use the report figures.")
}

#' Publication figure theme
#'
#' [ggpubr::theme_pubr()] sized for a figure drawn at its final printed width:
#' 7 pt text, 8 pt axis titles, heavier axis lines and ticks, and no title,
#' subtitle or caption. The font defaults to Arial. Where Arial is not
#' installed, fontconfig substitutes a metric-compatible font such as
#' Liberation Sans, so the layout is the same.
#'
#' @param base_size Text size in points. Default `7`, the smallest size Elsevier
#'   accepts in print and the largest Nature allows for non-label text.
#' @param base_family Font family. Default `"Arial"`.
#' @param legend Legend position. Default `"none"`.
#' @return A ggplot2 theme.
#' @export
ap_theme_pub <- function(base_size = 7, base_family = "Arial", legend = "none") {
  ap_pub_require("ggpubr")
  ggpubr::theme_pubr(base_size = base_size, base_family = base_family, legend = legend) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(size = base_size, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size + 1, colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = ggplot2::unit(1.2, "mm"),
      # Facets put their label on the left, outside the axis, where it reads as that
      # panel's y-axis title (see ap_pub_facet()).
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(size = base_size + 1, colour = "black",
                                         margin = ggplot2::margin(1, 1, 3, 1)),
      strip.text.y.left = ggplot2::element_text(size = base_size + 1, colour = "black",
                                                angle = 90, margin = ggplot2::margin(0, 2, 0, 0)),
      legend.text = ggplot2::element_text(size = base_size),
      legend.title = ggplot2::element_text(size = base_size + 1),
      legend.key.size = ggplot2::unit(3.5, "mm"),
      panel.spacing = ggplot2::unit(3, "mm"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 3, 2, 2, unit = "mm")
    )
}

#' Publication colour palette
#'
#' A [ggsci](https://CRAN.R-project.org/package=ggsci) journal palette. These
#' copy journal house styles and were not designed for colour vision
#' deficiency, unlike the Okabe-Ito palette the report figures use; check the
#' palette you choose against your group count.
#'
#' @param n Number of colours.
#' @param name ggsci palette: `"npg"` (default), `"lancet"`, `"nejm"`, `"jco"`
#'   or `"aaas"`.
#' @return A character vector of `n` hex colours.
#' @export
ap_pub_palette <- function(n, name = c("npg", "lancet", "nejm", "jco", "aaas")) {
  name <- match.arg(name)
  ap_pub_require("ggsci")
  pal <- getExportedValue("ggsci", paste0("pal_", name))()
  # Asking for more colours than a ggsci palette has returns NA with a warning.
  max_n <- sum(!is.na(suppressWarnings(pal(100))))
  ap_assert(n <= max_n,
            "The {name} palette has {max_n} colours; {n} groups need more. Choose another palette.")
  pal(n)
}

#' Options for publication figures
#'
#' The choices that differ between studies, collected once and passed to every
#' `ap_plot_*(publication = TRUE)` call. In the workflow they come from the
#' `publication:` block of the config. Everything else (column widths, font,
#' text sizes) is fixed by journal rules and is not an option.
#'
#' @param palette ggsci palette, see [ap_pub_palette()]. Default `"npg"`.
#' @param labels Named list mapping metadata variables to axis titles, e.g.
#'   `list(dx = "Diagnosis")`. A variable not listed is shown by its name.
#' @param capitalize_levels Capitalize the first letter of group levels on the
#'   axis ("adenoma" becomes "Adenoma"). Default `TRUE`. Does not touch a level
#'   renamed in `level_labels`, which is shown exactly as written.
#' @param level_labels Named list mapping group levels, as they appear in the
#'   metadata, to the names the figures and legends show, e.g.
#'   `list(cancer = "CRC", normal = "Healthy")`. A level not listed keeps its
#'   metadata value. Two levels may not share a name.
#' @return An `ap_pub_options` list.
#' @export
ap_pub_options <- function(palette = "npg", labels = list(), capitalize_levels = TRUE,
                           level_labels = list()) {
  palette <- match.arg(palette, c("npg", "lancet", "nejm", "jco", "aaas"))
  ap_assert(is.list(labels) && (length(labels) == 0L || !is.null(names(labels))),
            "`labels` must be a named list, e.g. list(dx = \"Diagnosis\").")
  ap_assert(isTRUE(capitalize_levels) || isFALSE(capitalize_levels),
            "`capitalize_levels` must be TRUE or FALSE.")
  ap_assert(is.list(level_labels) &&
              (length(level_labels) == 0L || (!is.null(names(level_labels)) &&
                                              all(nzchar(names(level_labels))))),
            "`level_labels` must be a named list, e.g. list(cancer = \"CRC\").")
  ap_assert(all(vapply(level_labels, function(v) is.character(v) && length(v) == 1L &&
                         nzchar(v), logical(1))),
            "Each entry of `level_labels` must be one non-empty name.")
  # Two levels renamed to one name would be merged into a single group without a word.
  shown <- unlist(level_labels)
  ap_assert(!anyDuplicated(shown),
            "`level_labels` gives two levels the same name ({paste(unique(shown[duplicated(shown)]), collapse = ', ')}); they would be drawn as one group.")
  structure(list(palette = palette, labels = labels, capitalize_levels = capitalize_levels,
                 level_labels = level_labels),
            class = "ap_pub_options")
}

#' @keywords internal
ap_pub_label <- function(var, pub) {
  lab <- pub$labels[[var]]
  if (is.null(lab)) var else as.character(lab)
}

# Group levels as figure labels: renamed where `level_labels` says so (exactly as written),
# otherwise capitalized when `capitalize_levels` is set.
#' @keywords internal
ap_pub_levels <- function(f, pub) {
  f <- as.factor(f)
  lv <- levels(f)
  mapped <- pub$level_labels[lv]
  renamed <- !vapply(mapped, is.null, logical(1))
  out <- lv
  if (isTRUE(pub$capitalize_levels)) out <- paste0(toupper(substr(lv, 1, 1)), substring(lv, 2))
  out[renamed] <- unlist(mapped[renamed])
  ap_assert(!anyDuplicated(out),
            "Two groups would be shown with the same name ({paste(unique(out[duplicated(out)]), collapse = ', ')}). Check `level_labels`.")
  levels(f) <- out
  f
}

# A group level as it reads inside a sentence (axis titles, headings, legends): the
# `level_labels` name when there is one, otherwise the metadata value unchanged.
#' @keywords internal
ap_pub_level_text <- function(v, pub) {
  v <- as.character(v)
  if (is.null(pub) || length(pub$level_labels) == 0L) return(v)
  hit <- v %in% names(pub$level_labels)
  v[hit] <- unlist(pub$level_labels[v[hit]])
  v
}

# Panel grid for a multi-panel figure: one row up to three panels, then two
# columns for four, three columns beyond. Each row is 60 mm tall. The size is
# attached to the plot so ap_save_figure() saves it at the size it was laid out for.
#' @keywords internal
ap_pub_layout <- function(n) {
  ncol <- if (n <= 3L) n else if (n == 4L) 2L else 3L
  nrow <- ceiling(n / ncol)
  list(ncol = ncol, nrow = nrow,
       width = if (n == 1L) "single" else "double",
       height = min(247, nrow * 60))
}

# facet_wrap with each panel's label as its y-axis title.
#' @keywords internal
ap_pub_facet <- function(facet, n) {
  lay <- ap_pub_layout(n)
  list(ggplot2::facet_wrap(stats::as.formula(paste("~", facet)), ncol = lay$ncol,
                           scales = "free_y", strip.position = "left"),
       lay)
}
