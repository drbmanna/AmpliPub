# Publication figures: journal sizes, three formats, no statistics on the panel,
# and the statistics and caveats carried into the legend file instead.

skip_pub <- function() {
  for (p in c("ggpubr", "ggsci", "svglite")) skip_if_not_installed(p)
}

has_geom <- function(p, cls) any(vapply(p$layers, function(l) inherits(l$geom, cls), logical(1)))

# Width read back from the files. The PDF is not read: cairo writes its page
# dictionary into a compressed object stream, and no PDF parser is a dependency.
svg_width_pt <- function(path) {
  head <- paste(readLines(path, n = 5L, warn = FALSE), collapse = " ")
  as.numeric(sub(".*<svg[^>]* width='([0-9.]+)pt'.*", "\\1", head))
}
png_width_px <- function(path) {
  b <- readBin(path, "raw", 24L)
  sum(as.integer(b[17:20]) * 256^(3:0))
}

test_that("the publication alpha figure has no caption, no statistics and the pub theme", {
  skip_pub()
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = c("q0", "q1"), n_iter = 2L)
  p <- ap_plot_alpha(a, "group", publication = TRUE)
  expect_null(p$labels$caption)
  expect_null(p$labels$subtitle)
  expect_false(has_geom(p, "GeomText"))
  expect_true(inherits(p$theme$plot.caption, "element_blank"))
  expect_equal(p$theme$axis.text$size, 7)
  expect_silent(ggplot2::ggplot_build(p))
  # The report figure is unchanged: it still carries the statistics.
  expect_true(has_geom(ap_plot_alpha(a, "group"), "GeomText"))
})

test_that("ap_save_figure writes PDF, SVG, PNG and legend at the journal width", {
  skip_pub()
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = "q0", n_iter = 2L)
  p <- ap_plot_alpha(a, "group", publication = TRUE)
  dir <- withr::local_tempdir()
  files <- ap_save_figure(p, "alpha", dir, width = "single", height = 60,
                          legend = ap_alpha_legend(a))
  expect_setequal(basename(files), c("alpha.pdf", "alpha.svg", "alpha.png", "alpha_legend.txt"))
  expect_true(all(file.size(files) > 0))
  expect_equal(svg_width_pt(file.path(dir, "alpha.svg")), 89 / 25.4 * 72, tolerance = 0.01)
  expect_equal(png_width_px(file.path(dir, "alpha.png")), round(89 / 25.4 * 600), tolerance = 1)
  expect_identical(readBin(file.path(dir, "alpha.pdf"), "raw", 5L), charToRaw("%PDF-"))
  svg <- readLines(file.path(dir, "alpha.svg"))
  expect_true(any(grepl("font-family: Arial, \"Liberation Sans\", sans-serif;", svg, fixed = TRUE)))
  expect_false(any(grepl("font-family: \"Liberation Sans\";", svg, fixed = TRUE)))
})

test_that("the panel grid follows the number of panels and sets the saved size", {
  lay <- lapply(1:9, ap_pub_layout)
  expect_equal(vapply(lay, `[[`, numeric(1), "ncol"), c(1, 2, 3, 2, 3, 3, 3, 3, 3))
  expect_equal(vapply(lay, `[[`, numeric(1), "nrow"), c(1, 1, 1, 2, 2, 2, 3, 3, 3))
  expect_equal(lay[[1]]$width, "single")
  expect_equal(lay[[5]]$height, 120)
  expect_true(all(vapply(lapply(1:30, ap_pub_layout), `[[`, numeric(1), "height") <= 247))

  skip_pub()
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  p <- ap_plot_alpha(ap_alpha(x, metrics = c("q0", "q1", "q2", "evenness"), n_iter = 2L),
                     "group", publication = TRUE)
  expect_equal(attr(p, "ap_pub_size"), list(width = "double", height = 120))
})

test_that("publication options label the variable, capitalize levels and pick the palette", {
  skip_pub()
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = "q0", n_iter = 2L)
  pub <- ap_pub_options(palette = "lancet", labels = list(group = "Treatment arm"))
  p <- ap_plot_alpha(a, "group", publication = TRUE, pub = pub)
  expect_equal(p$labels$x, "Treatment arm")
  b <- ggplot2::ggplot_build(p)
  lv <- levels(b$plot$data$grp)
  expect_true(all(substr(lv, 1, 1) == toupper(substr(lv, 1, 1))))
  expect_setequal(unique(b$data[[2]]$fill), ap_pub_palette(length(lv), "lancet"))

  p0 <- ap_plot_alpha(a, "group", publication = TRUE,
                      pub = ap_pub_options(capitalize_levels = FALSE))
  expect_equal(p0$labels$x, "group")
  expect_equal(levels(p0$data$grp), levels(factor(a$metadata$group)))

  expect_error(ap_pub_options(palette = "rainbow"), "should be one of")
  expect_error(ap_pub_options(labels = list("x")), "named list")
})

test_that("ap_save_figure refuses sizes a journal page cannot hold", {
  skip_pub()
  p <- ggplot2::ggplot()
  dir <- withr::local_tempdir()
  expect_error(ap_save_figure(p, "x", dir, width = "triple", height = 50), "must be one of")
  expect_error(ap_save_figure(p, "x", dir, width = 200, height = 50), "between 0 and 183")
  expect_error(ap_save_figure(p, "x", dir, width = "single", height = 300), "247")
})

test_that("ap_pub_palette refuses more groups than the palette has", {
  skip_pub()
  expect_length(ap_pub_palette(3, "npg"), 3)
  expect_error(ap_pub_palette(50, "npg"), "has 10 colours")
})

test_that("a table already at one depth is not called unrarefied", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  counts <- SummarizedExperiment::assay(x, "counts")
  set.seed(1)
  SummarizedExperiment::assay(x, "counts") <- ap_rarefy_matrix(counts, 500)
  a <- ap_alpha(x, metrics = "q0", rarefy = FALSE)
  expect_equal(a$common_depth, 500)
  expect_match(ap_alpha_legend(a)[1], "common depth of 500 reads")
  expect_no_match(ap_alpha_legend(a)[1], "Not rarefied")

  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  expect_equal(b$common_depth, 500)
  expect_no_match(paste(cli::cli_fmt(print(b)), collapse = "\n"), "Not rarefied")
})

test_that("a table at unequal depths keeps the not-rarefied caveat", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  # The fixture has every library at 5,000 reads; double one to make depths unequal.
  counts <- SummarizedExperiment::assay(x, "counts")
  counts[, 1] <- counts[, 1] * 2L
  SummarizedExperiment::assay(x, "counts") <- counts
  a <- suppressWarnings(ap_alpha(x, metrics = "q0", rarefy = FALSE))
  expect_true(is.na(a$common_depth))
  expect_match(ap_alpha_legend(a)[1], "Not rarefied")
})

test_that("the legend carries each metric's test, n, q and effect size", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  a <- ap_alpha(x, metrics = c("q0", "q1"), n_iter = 2L)
  t <- ap_alpha_test(a, "group", seed = 1L)
  leg <- ap_alpha_legend(a, t)
  expect_length(leg, 3)
  expect_match(leg[2], "^Hill q0 \\(observed richness\\): .+, n = \\d+, q = .+, .+ = -?[0-9.]+")
})

ord_fixture <- function() {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  list(ord = ap_ordinate(b, "bray_curtis", method = "pcoa"),
       pn = ap_permanova(b, "group", permutations = 99L, seed = 1L))
}

test_that("the publication ordination drops the verdict but keeps the numbers", {
  skip_pub()
  f <- ord_fixture()
  p <- ap_plot_ordination(f$ord, group = "group", permanova = f$pn, publication = TRUE,
                          pub = ap_pub_options(labels = list(group = "Arm")))
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  grobs <- lapply(Filter(function(l) inherits(l$geom, "GeomCustomAnn"), p$layers),
                  function(l) l$geom_params$grob)
  txt <- vapply(grobs, function(g) paste(deparse(g$label), collapse = ""), character(1))
  expect_length(txt, 2)
  expect_match(txt[1], "italic(R)^2", fixed = TRUE)
  expect_match(txt[1], "italic(F)", fixed = TRUE)
  expect_match(txt[2], "betadisper", fixed = TRUE)
  # Baselines 2.8 mm apart, whatever each line's height.
  ys <- vapply(grobs, function(g) grid::convertY(grid::unit(1, "npc") - g$y, "mm", valueOnly = TRUE),
               numeric(1))
  expect_equal(unname(diff(ys)), 2.8)
  verdicts <- f$pn$interpretation$verdict
  expect_false(any(vapply(verdicts, function(v) any(grepl(v, txt, ignore.case = TRUE)), logical(1))))
  # Every line is a plotmath expression, so the PDF gets a real superscript.
  expect_true(all(vapply(grobs, function(g) is.call(g$label), logical(1))))
  # Ellipses: one filled polygon layer at low opacity, one outline.
  polys <- Filter(function(l) inherits(l$geom, "GeomPolygon"), p$layers)
  expect_length(polys, 1)
  expect_lt(polys[[1]]$aes_params$alpha, 0.2)
  expect_equal(p$scales$get_scales("colour")$name, "Arm")
  expect_equal(p$scales$get_scales("fill")$name, "Arm")
  size <- attr(p, "ap_pub_size")
  expect_equal(size$width, "single")
  expect_true(size$height >= 50 && size$height <= 247)
  expect_equal(unname(vapply(grobs, function(g) g$gp$fontsize, numeric(1))), c(6, 6))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("the ordination legend carries the verdict and its reasoning", {
  f <- ord_fixture()
  leg <- ap_ordination_legend(f$ord, f$pn, "group")
  expect_match(leg[1], "^PCoA on bray_curtis")
  expect_true(any(startsWith(leg, "PERMANOVA (99 permutations): R2 = ")))
  expect_true(any(grepl(paste0("^Verdict: ", f$pn$interpretation$verdict[1]), leg)))
})

test_that("the panel statistics parse, including a missing dispersion test", {
  pn <- list(results = data.frame(metric = "m", term = "t", R2 = 0.1, pseudo_F = 2, p = 0.001),
             dispersion = data.frame(metric = "m", term = "t", dispersion_p = NA_real_))
  lines <- ap_permanova_plotmath(pn, "m", "t")
  expect_match(lines[1], "italic(p) == '0.001'", fixed = TRUE)
  expect_equal(lines[2], "'betadisper:' ~ '-'")
  for (t in lines) expect_silent(parse(text = t))
  expect_null(ap_permanova_plotmath(pn, "m", "other"))
})
