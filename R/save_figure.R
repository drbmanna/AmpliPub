# Saving publication figures.
#
# Widths are the column widths both Nature and Elsevier accept (checked
# 2026-09-18): Nature 89 / 120 / 183 mm, Elsevier 90 / 140 / 190 mm. The
# smaller of each pair fits both. Nature's full page depth, 247 mm, caps the
# height.

ap_figure_widths <- c(single = 89, onehalf = 120, double = 183)

#' Save a publication figure
#'
#' Writes the figure as PDF (vector, fonts embedded), SVG (vector, text left
#' editable) and PNG, at a journal column width, plus `<name>_legend.txt`
#' holding the statistics and caveats that the publication figure leaves off
#' the panel. A figure is never exported without the caveats its statistic came
#' with; they move from the panel to the legend.
#'
#' In the SVG the font is written as `Arial, "Liberation Sans", sans-serif`, so
#' it shows in Arial wherever Arial is installed. The PDF embeds whatever font
#' drew it, which on a machine without Arial is its metric-compatible
#' substitute.
#'
#' @param plot A ggplot object.
#' @param name File name without extension.
#' @param dir Output directory, created if missing.
#' @param width `"single"` (89 mm), `"onehalf"` (120 mm), `"double"` (183 mm),
#'   or a number of millimetres. Defaults to the size the `ap_plot_*` function
#'   laid the figure out for.
#' @param height Height in millimetres, at most 247. Default as for `width`.
#' @param legend Character vector for the legend file. `NULL` writes none.
#' @param dpi PNG resolution. Default `600`.
#' @return The paths written, invisibly.
#' @export
ap_save_figure <- function(plot, name, dir, width = NULL, height = NULL, legend = NULL,
                           dpi = 600) {
  ap_assert(inherits(plot, "ggplot"), "`plot` must be a ggplot, not {class(plot)[1]}.")
  size <- attr(plot, "ap_pub_size")
  width <- width %||% size$width
  height <- height %||% size$height
  ap_assert(!is.null(width) && !is.null(height),
            "Give `width` and `height`; this plot does not carry a publication size.")
  if (is.character(width)) {
    ap_assert(length(width) == 1L && width %in% names(ap_figure_widths),
              "`width` must be one of {paste(names(ap_figure_widths), collapse = ', ')} or a number of mm.")
    width <- ap_figure_widths[[width]]
  }
  ap_assert(is.numeric(width) && width > 0 && width <= max(ap_figure_widths),
            "`width` must be between 0 and {max(ap_figure_widths)} mm.")
  ap_assert(is.numeric(height) && height > 0 && height <= 247,
            "`height` must be between 0 and 247 mm, the full page depth.")
  ap_pub_require("svglite")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  path <- function(ext) file.path(dir, paste0(name, ext))
  save <- function(file, device, ...) {
    ggplot2::ggsave(file, plot, width = width, height = height, units = "mm",
                    device = device, ...)
  }
  save(path(".pdf"), grDevices::cairo_pdf)
  save(path(".svg"), svglite::svglite)
  save(path(".png"), if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png",
       dpi = dpi, bg = "white")

  svg <- readLines(path(".svg"), warn = FALSE)
  svg <- gsub("font-family: [^;\"]*(\"[^\"]*\")?[^;]*;", "font-family: Arial, \"Liberation Sans\", sans-serif;", svg)
  writeLines(svg, path(".svg"))

  files <- path(c(".pdf", ".svg", ".png"))
  if (!is.null(legend)) {
    writeLines(legend, path("_legend.txt"))
    files <- c(files, path("_legend.txt"))
  }
  bad <- files[!file.exists(files) | file.size(files) == 0]
  ap_assert(length(bad) == 0L, "Figure files missing or empty: {paste(basename(bad), collapse = ', ')}.")
  invisible(files)
}
