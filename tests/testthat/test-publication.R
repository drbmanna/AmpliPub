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
  # Caption, the two naming notes (Hill numbers; why no Chao1), then one line per metric.
  expect_length(leg, 5)
  expect_match(leg[2], "^Shannon and inverse Simpson are the Hill numbers of order 1 and 2")
  expect_match(leg[3], "^Chao1 and ACE are not reported")
  expect_match(leg[4], "^Richness: .+, n = \\d+, q = .+, .+ = -?[0-9.]+")
  expect_match(leg[5], "^Shannon: ")
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

test_that("the publication dispersion carries both tests on the panel, not the verdict", {
  skip_pub()
  f <- ord_fixture()
  p <- ap_plot_dispersion(f$pn, metric = "bray_curtis", term = "group", publication = TRUE,
                          pub = ap_pub_options(labels = list(group = "Arm")))
  expect_null(p$labels$subtitle)
  expect_null(p$labels$caption)
  grobs <- lapply(Filter(function(l) inherits(l$geom, "GeomCustomAnn"), p$layers),
                  function(l) l$geom_params$grob)
  txt <- vapply(grobs, function(g) paste(deparse(g$label), collapse = ""), character(1))
  # PERMANOVA on top, then betadisper, then the spread ratio.
  expect_length(txt, 3)
  expect_match(txt[1], "PERMANOVA", fixed = TRUE)
  expect_match(txt[1], "italic(R)^2", fixed = TRUE)
  expect_match(txt[2], "betadisper", fixed = TRUE)
  expect_match(txt[3], "spread ratio", fixed = TRUE)
  # Baselines 2.8 mm apart, as on the publication ordination.
  ys <- vapply(grobs, function(g) grid::convertY(grid::unit(1, "npc") - g$y, "mm", valueOnly = TRUE),
               numeric(1))
  expect_equal(unname(diff(ys)), c(2.8, 2.8))
  verdicts <- f$pn$interpretation$verdict
  expect_false(any(vapply(verdicts, function(v) any(grepl(v, txt, ignore.case = TRUE)), logical(1))))
  expect_true(all(vapply(grobs, function(g) is.call(g$label), logical(1))))
  expect_equal(unname(vapply(grobs, function(g) g$gp$fontsize, numeric(1))), c(6, 6, 6))
  expect_equal(as.character(p$labels$x), "Arm")
  size <- attr(p, "ap_pub_size")
  expect_equal(size$width, "single")
  expect_true(size$height >= 50 && size$height <= 247)
  # The report figure is unchanged: it still carries a subtitle.
  expect_false(is.null(ap_plot_dispersion(f$pn, "bray_curtis", "group")$labels$subtitle))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("the dispersion legend carries the method note and the verdict", {
  f <- ord_fixture()
  leg <- ap_dispersion_legend(f$pn, "bray_curtis", "group")
  expect_match(leg[1], "^Distance from each sample to its group centroid on bray_curtis")
  expect_true(any(grepl("PERMANOVA assumes", leg)))
  expect_true(any(startsWith(leg, "betadisper: F = ")))
  expect_true(any(grepl(paste0("^Verdict: ", f$pn$interpretation$verdict[1]), leg)))
})

test_that("the dispersion panel statistics parse, and are dropped when betadisper is missing", {
  pn <- list(results = data.frame(metric = "m", term = "t", R2 = 0.006, pseudo_F = 1.5, p = 0.005),
             dispersion = data.frame(metric = "m", term = "t", dispersion_F = 3.6,
                                     dispersion_p = 0.034, max_centroid_ratio = 1.04))
  lines <- ap_dispersion_plotmath(pn, "m", "t")
  expect_length(lines, 3)
  expect_match(lines[1], "PERMANOVA", fixed = TRUE)
  expect_match(lines[3], "1.04", fixed = TRUE)
  for (t in lines) expect_silent(parse(text = t))
  pn$dispersion$dispersion_p <- NA_real_
  expect_null(ap_dispersion_plotmath(pn, "m", "t"))
})

test_that("publication taxon labels drop the rank prefix and keep every suffix", {
  lab <- ap_pub_taxon_labels(c("g__Blautia_A_141781", "g__Faecalibacterium",
                               "f__Lachnospiraceae (unassigned genus)", "Other"))
  txt <- vapply(lab, function(e) paste(deparse(e), collapse = ""), character(1))
  expect_equal(unname(txt[1]), "italic(\"Blautia_A_141781\")")
  expect_equal(unname(txt[2]), "italic(\"Faecalibacterium\")")
  expect_equal(unname(txt[3]), "italic(\"Lachnospiraceae\") ~ \"(unassigned genus)\"")
  expect_equal(unname(txt[4]), "\"Other\"")
  expect_equal(names(lab), c("g__Blautia_A_141781", "g__Faecalibacterium",
                             "f__Lachnospiraceae (unassigned genus)", "Other"))
})

test_that("two lineages differing only by a numeric suffix keep distinct labels", {
  lab <- ap_pub_taxon_labels(c("g__Blautia_A_141780", "g__Blautia_A_141781"))
  txt <- vapply(lab, function(e) paste(deparse(e), collapse = ""), character(1))
  expect_equal(anyDuplicated(txt), 0L)
})

test_that("the label guard fires when removing prefixes would merge two taxa", {
  expect_error(ap_pub_taxon_labels(c("g__Bacteroides", "f__Bacteroides")),
               "same label")
})

test_that("the taxa palette has no greys and grows past twenty colours", {
  skip_if_not_installed("ggsci")
  pal <- ap_pub_taxa_palette(30)
  expect_length(pal, 30L)
  expect_equal(anyDuplicated(toupper(pal)), 0L)
  rgb <- grDevices::col2rgb(pal)
  expect_true(all(apply(rgb, 2, function(z) diff(range(z))) >= 16))
})

test_that("the publication composition bar shows group means with italic taxa", {
  skip_pub()
  x <- ap_fixture_object(tree = FALSE)
  p <- ap_plot_taxa_bar(x, rank = "genus", n = 3L, group = "group", publication = TRUE)
  expect_null(p$labels$caption)
  expect_equal(attr(p, "ap_pub_size")$width, "onehalf")
  b <- ggplot2::ggplot_build(p)
  # One stacked bar per group, each summing to 1.
  sums <- tapply(b$data[[1]]$y - b$data[[1]]$ymin, b$data[[1]]$x, sum)
  expect_equal(unname(as.vector(sums)), rep(1, length(sums)), tolerance = 1e-8)
  expect_error(ap_plot_taxa_bar(x, rank = "genus", n = 3L, publication = TRUE),
               "needs `group`")
})

test_that("the publication heatmap has one column per group and no Other row", {
  skip_pub()
  x <- ap_fixture_object(tree = FALSE)
  h <- ap_plot_taxa_heatmap(x, rank = "genus", n = 3L, group = "group", publication = TRUE)
  d <- ggplot2::ggplot_build(h)$data[[1]]
  n_grp <- length(unique(SummarizedExperiment::colData(x)$group))
  expect_equal(length(unique(d$x)), n_grp)
  expect_false("Other" %in% levels(h$data$taxon))
  expect_equal(nrow(d), n_grp * nlevels(h$data$taxon))
  expect_equal(h$labels$y, "Genus")
})

test_that("without a reference each heatmap row is centred on the group average, not scaled", {
  x <- ap_fixture_object(tree = FALSE)
  df <- ap_top_taxa(x, n = 3L, group = "group", rank = "genus")
  m <- ap_taxa_clr_contrast(x, df, "group", "genus", NULL, 0.5)
  expect_equal(unname(rowSums(m)), rep(0, nrow(m)), tolerance = 1e-10)
  # Centred only: the spread of a row is the spread of its group means, not 1.
  clr <- ap_normalize(ap_collapse(x, "genus"), method = "clr", pseudocount = 0.5)
  g <- SummarizedExperiment::colData(x)$group
  t1 <- rownames(m)[1]
  raw <- tapply(clr[t1, ], g[match(colnames(clr), colnames(x))], mean)
  expect_equal(unname(m[t1, names(raw)]), as.vector(raw - mean(raw)), tolerance = 1e-10)
})

test_that("with a reference the heatmap shows differences from it and drops its column", {
  x <- ap_fixture_object(tree = FALSE)
  df <- ap_top_taxa(x, n = 3L, group = "group", rank = "genus")
  lv <- sort(unique(as.character(SummarizedExperiment::colData(x)$group)))
  m0 <- ap_taxa_clr_contrast(x, df, "group", "genus", NULL, 0.5)
  m1 <- ap_taxa_clr_contrast(x, df, "group", "genus", lv[1], 0.5)
  expect_false(lv[1] %in% colnames(m1))
  expect_equal(m1[, lv[2]], m0[, lv[2]] - m0[, lv[1]], tolerance = 1e-10)
  expect_error(ap_taxa_clr_contrast(x, df, "group", "genus", "no_such_level", 0.5),
               "must be a level")
})

test_that("the composition legend states selection, Other's share and group sizes", {
  x <- ap_fixture_object(tree = FALSE)
  df <- ap_top_taxa(x, n = 2L, group = "group", rank = "genus")
  leg <- ap_taxa_legend(x, "genus", "group", n = 2L, type = "bar")
  expect_match(leg[2], sprintf("(%d taxa)", attr(df, "n_kept")), fixed = TRUE)
  if (attr(df, "n_pooled") > 0L) {
    expect_match(leg[3], sprintf("%.1f%%", 100 * attr(df, "pooled_mean_abundance")), fixed = TRUE)
  }
  expect_match(leg[5], "^n = ")
  expect_match(ap_taxa_legend(x, "genus", "group", n = 2L, type = "heatmap")[1], "average of the group means")
})

# Differential abundance publication figures ---------------------------------------------------

fake_da <- function() {
  f <- c("aaaaaa11", "bbbbbb22", "cccccc33", "dddddd44")
  lab <- c("g__Prevotella", "g__Prevotella", "g__Parvimonas", "f__Lachnospiraceae (unassigned genus)")
  rows <- expand.grid(i = 1:4, method = c("ancombc2", "linda"), stringsAsFactors = FALSE)
  r <- data.frame(feature = f[rows$i], method = rows$method, contrast = "b_vs_a",
                  effect = c(1.2, -0.8, 2, 0.1, 1.1, -0.1, 1.5, 0.05),
                  effect_scale = ifelse(rows$method == "ancombc2", "log fold change (natural log)",
                                        "log2 fold change"),
                  se = 0.2, statistic = 5, p = 1e-4,
                  p_adj = c(1e-4, 1e-3, 1e-8, 0.6, 1e-3, 0.5, 1e-5, 0.9),
                  note = NA_character_, stringsAsFactors = FALSE)
  r$passed_sensitivity <- ifelse(r$method == "ancombc2", c(FALSE, TRUE, TRUE, TRUE), NA)
  below <- r$p_adj < 0.05
  r$failed_sensitivity <- below & r$passed_sensitivity %in% FALSE
  r$significant <- below & !r$failed_sensitivity
  r$prevalence <- 0.5
  r$taxon_label <- lab[rows$i]
  structure(list(results = r, group = "group", reference = "a", levels = c("a", "b"),
                 methods = c("ancombc2", "linda"), skipped = character(0), prv_cut = 0.1,
                 alpha = 0.05, p_adj_method = "BH", n_features = 4L, n_samples = 20L),
            class = "ap_da")
}

test_that("ASVs sharing a genus get a short ID, others keep the plain label", {
  lab <- ap_da_feature_labels(c("aaaaaa11", "bbbbbb22", "cccccc33"),
                              c("g__Prevotella", "g__Prevotella", "g__Parvimonas"))
  expect_equal(lab, c("g__Prevotella (ASV aaaaaa)", "g__Prevotella (ASV bbbbbb)", "g__Parvimonas"))
  expect_error(ap_da_feature_labels(c("aaaaaa11", "aaaaaa22"), c("g__X", "g__X")), "not unique")
})

test_that("publication feature labels italicise the name and keep notes upright", {
  e <- ap_pub_feature_expr(c("g__Prevotella (ASV aaaaaa)", "f__Lachnospiraceae (unassigned genus)",
                             "g__Parvimonas"))
  txt <- vapply(e, function(x) paste(deparse(x), collapse = ""), character(1))
  expect_equal(unname(txt), c("italic(\"Prevotella\") ~ \"(ASV aaaaaa)\"",
                              "italic(\"Lachnospiraceae\") ~ \"(unassigned genus)\"",
                              "italic(\"Parvimonas\")"))
})

test_that("the publication concordance has a row per called feature, fragile calls as open circles", {
  skip_pub()
  da <- fake_da()
  cc <- ap_da_concordance(da)
  p <- ap_plot_concordance(cc, publication = TRUE)
  # Called by some method: aaaa (linda), bbbb (ancombc2), cccc (both). dddd by none.
  expect_setequal(unique(as.character(p$data$feature)), c("aaaaaa11", "bbbbbb22", "cccccc33"))
  fail <- p$data[p$data$kind == "Failed sensitivity analysis", ]
  expect_equal(fail$feature, "aaaaaa11")
  # A structural zero is drawn as a triangle coloured by its direction.
  da2 <- fake_da()
  sz <- da2$results[1, ]
  sz$feature <- "eeeeee55"; sz$taxon_label <- "g__Fusobacterium"; sz$effect <- NA
  sz$p_adj <- NA; sz$p <- NA; sz$failed_sensitivity <- FALSE; sz$passed_sensitivity <- NA
  sz$significant <- TRUE
  da2$results$structural_zero <- FALSE
  da2$results$direction <- sign(da2$results$effect)
  sz$structural_zero <- TRUE; sz$direction <- 1
  da2$results <- rbind(da2$results, sz)
  p2 <- ap_plot_concordance(ap_da_concordance(da2), publication = TRUE)
  tri <- p2$data[p2$data$kind == "Absent from one group", ]
  expect_equal(tri$feature, "eeeeee55")
  expect_equal(tri$colour, "Enriched")
  expect_silent(ggplot2::ggplot_build(p2))
  expect_silent(ggplot2::ggplot_build(p))
  expect_equal(attr(p, "ap_pub_size")$width, "single")
  # Only the exceptions are keyed, once each: a filled dot needs no key, and an absent
  # triangle layer with no rows must not add a second circle to the failed key.
  keys <- ggplot2::get_guide_data(p, "shape")
  expect_equal(keys$.label, "Failed sensitivity analysis")
  expect_setequal(ggplot2::get_guide_data(p2, "shape")$.label,
                  c("Absent from one group", "Failed sensitivity analysis"))
})

test_that("the publication volcano titles each panel with its scale, without nested brackets", {
  skip_pub()
  v <- ap_plot_volcano(fake_da(), publication = TRUE)
  expect_true(all(grepl("^(ANCOM-BC2|LinDA)\n[(]", levels(v$data$panel))))
  expect_false(any(grepl("[(].*[(]", levels(v$data$panel))))
  expect_equal(attr(v, "ap_pub_size")$width, "double")
  expect_silent(ggplot2::ggplot_build(v))
})

test_that("the DA legend states counts, the failed-sensitivity number and full IDs", {
  cc <- ap_da_concordance(fake_da())
  leg <- ap_da_legend(cc, "concordance")
  expect_true(any(grepl("Called: ANCOM-BC2 2, LinDA 2.", leg, fixed = TRUE)))
  expect_true(any(grepl("^1 ANCOM-BC2 call had adjusted p < 0.05 but failed", leg)))
  expect_true(any(grepl("g__Prevotella (ASV aaaaaa) = aaaaaa11", leg, fixed = TRUE)))
})

# Primary-method figure ---------------------------------------------------------------------------

# n enriched and n depleted ANCOM-BC2 calls with distinct effects, for the per-direction cap.
many_calls_da <- function(n = 30L) {
  f <- sprintf("f%03d%s", seq_len(2L * n), strrep("0", 30))
  eff <- c(seq(0.5, 3, length.out = n), -seq(0.5, 3, length.out = n))
  r <- data.frame(feature = f, method = "ancombc2", contrast = "b_vs_a", effect = eff,
                  effect_scale = "log fold change (natural log)", se = 0.1, statistic = 5,
                  p = 1e-4, p_adj = 1e-3, note = NA_character_, passed_sensitivity = TRUE,
                  failed_sensitivity = FALSE, significant = TRUE, prevalence = 0.5,
                  taxon_label = sprintf("g__G%03d", seq_len(2L * n)), stringsAsFactors = FALSE)
  structure(list(results = r, group = "group", reference = "a", levels = c("a", "b"),
                 methods = "ancombc2", skipped = character(0), prv_cut = 0.1, alpha = 0.05,
                 p_adj_method = "BH", n_features = 2L * n, n_samples = 20L), class = "ap_da")
}

test_that("the primary figure draws only the primary method's robust calls, enriched first", {
  skip_pub()
  cc <- ap_da_concordance(fake_da())
  p <- ap_plot_da_primary(cc, "ancombc2")
  # ancombc2 called cccc (+2) and bbbb (-0.8); aaaa failed its sensitivity analysis.
  expect_equal(p$data$feature, c("cccccc33", "bbbbbb22"))
  expect_true(p$data$y[p$data$feature == "cccccc33"] > p$data$y[p$data$feature == "bbbbbb22"])
  heads <- attr(p, "ap_heads")
  expect_equal(heads$text, c("Enriched in b", "Depleted in b"))
  # Enriched heading right-aligned at the right edge, depleted left-aligned at the left.
  expect_equal(heads$x, c(Inf, -Inf))
  expect_equal(heads$hjust, c(1, 0))
  # Headings sit on their own rows, above each block.
  expect_equal(heads$y[1], max(p$data$y[p$data$direction > 0]) + 1)
  expect_equal(heads$y[2], max(p$data$y[p$data$direction < 0]) + 1)
  expect_silent(ggplot2::ggplot_build(p))
  expect_equal(attr(p, "ap_pub_size")$width, "single")
})

test_that("the primary figure's interval is estimate +/- 1.96 SE", {
  skip_pub()
  p <- ap_plot_da_primary(ap_da_concordance(fake_da()), "ancombc2")
  bars <- ggplot2::layer_data(p, 3L)
  z <- stats::qnorm(0.975)
  expect_equal(sort(bars$xmin), sort(c(2, -0.8) - z * 0.2))
  expect_equal(sort(bars$xmax), sort(c(2, -0.8) + z * 0.2))
})

test_that("the primary figure keeps the 25 largest effects in each direction", {
  skip_pub()
  cc <- ap_da_concordance(many_calls_da(30L))
  p <- ap_plot_da_primary(cc, "ancombc2")
  expect_equal(sum(p$data$direction > 0), 25L)
  expect_equal(sum(p$data$direction < 0), 25L)
  # The five smallest effects in each direction (0.5 upward) are the ones left out.
  expect_equal(min(abs(p$data$effect)), sort(seq(0.5, 3, length.out = 30))[6])
  expect_equal(unname(attr(p, "ap_n_called")), c(30L, 30L))
  leg <- ap_da_legend(cc, "primary", "ancombc2")
  expect_true(any(grepl("the 25 largest of 30 enriched and 25 largest of 30 depleted calls in b",
                        leg, fixed = TRUE)))
  expect_equal(nrow(ap_plot_da_primary(cc, "ancombc2", max_features = 3L)$data), 6L)
})

test_that("the primary figure refuses what it cannot draw honestly", {
  skip_pub()
  cc <- ap_da_concordance(fake_da())
  expect_error(ap_plot_da_primary(cc, "maaslin2"), "one of the methods that ran")
  # A primary method that called nothing has no figure.
  da <- fake_da()
  da$results$significant[da$results$method == "linda"] <- FALSE
  expect_error(ap_plot_da_primary(ap_da_concordance(da), "linda"), "called no feature")
  # ALDEx2's stored se is diff.win, a dispersion; no interval can be built from it.
  da2 <- fake_da()
  da2$results$method[da2$results$method == "linda"] <- "aldex2"
  da2$methods <- c("ancombc2", "aldex2")
  expect_error(ap_plot_da_primary(ap_da_concordance(da2), "aldex2"), "diff.win")
  # The abundance view is for a few chosen calls, and only calls.
  expect_error(ap_plot_da_primary(cc, type = "abundance"), "1 to 8 feature IDs")
  expect_error(ap_plot_da_primary(cc, features = "aaaaaa11"), "not called by ANCOM-BC2")
})

test_that("the primary legend states the method, interval, robust rule and what is shown", {
  cc <- ap_da_concordance(fake_da())
  leg <- ap_da_legend(cc, "primary", "ancombc2")
  expect_true(any(grepl("^Primary method ANCOM-BC2[.] Points: log fold change [(]natural log[)]", leg)))
  expect_true(any(grepl("1.96 SE, not adjusted for multiple testing", leg, fixed = TRUE)))
  expect_true(any(grepl("1 with adjusted p < 0.05 did not and is not shown.", leg, fixed = TRUE)))
  expect_true(any(grepl("All 2 calls shown (1 enriched, 1 depleted in b).", leg, fixed = TRUE)))
  expect_false(any(grepl("g__", leg)))
})

test_that("the abundance view uses the object DA ran on and counts detections per group", {
  skip_pub()
  x <- ap_fixture_object()
  counts <- SummarizedExperiment::assay(x, "counts")
  r <- data.frame(feature = c("ASV01", "ASV02"), method = "ancombc2", contrast = "b_vs_a",
                  effect = c(2, 1.5), effect_scale = "log fold change (natural log)", se = 0.2,
                  statistic = 5, p = 1e-4, p_adj = 1e-3, note = NA_character_,
                  passed_sensitivity = TRUE, failed_sensitivity = FALSE, significant = TRUE,
                  prevalence = 1, taxon_label = c("g__Genus01", "g__Genus02"),
                  stringsAsFactors = FALSE)
  da <- structure(list(results = r, group = "group", reference = "a", levels = c("a", "b"),
                       methods = "ancombc2", skipped = character(0), prv_cut = 0.1,
                       alpha = 0.05, p_adj_method = "BH", n_features = nrow(counts),
                       n_samples = ncol(counts)), class = "ap_da")
  cc <- ap_da_concordance(da)
  p <- ap_plot_da_primary(cc, type = "abundance", x = x, features = c("ASV01", "ASV02"))
  expect_silent(ggplot2::ggplot_build(p))
  # Relative abundance in percent, from all features' counts in each sample.
  v <- p$data$value[p$data$feature == "ASV01"]
  expect_equal(sort(v), sort(100 * counts["ASV01", ] / colSums(counts)), ignore_attr = TRUE)
  det <- ggplot2::layer_data(p, 3L)$label
  expect_true(all(grepl("^[0-9]+/12$", det)))
  # The wrong object is refused rather than drawn.
  expect_error(ap_plot_da_primary(cc, type = "abundance", x = x[, 1:10], features = "ASV01"),
               "Pass the object")
})

# Level labels ---------------------------------------------------------------------------------

test_that("level_labels rename levels exactly, and unlisted levels are capitalized", {
  pub <- ap_pub_options(level_labels = list(cancer = "CRC", normal = "healthy"))
  f <- ap_pub_levels(factor(c("normal", "cancer", "adenoma")), pub)
  # Renamed levels are shown as written ("healthy" stays lower case); the rest capitalized.
  expect_equal(levels(f), c("Adenoma", "CRC", "healthy"))
  expect_equal(as.character(f), c("healthy", "CRC", "Adenoma"))
  expect_equal(ap_pub_level_text(c("cancer", "adenoma"), pub), c("CRC", "adenoma"))
  expect_equal(ap_pub_level_text("cancer", NULL), "cancer")
  expect_error(ap_pub_options(level_labels = list("CRC")), "named list")
  expect_error(ap_pub_options(level_labels = list(cancer = "")), "non-empty")
  # Two levels renamed to one name would silently merge two groups.
  expect_error(ap_pub_options(level_labels = list(cancer = "Case", adenoma = "Case")),
               "same name")
  # A label that collides with another level's capitalized name is caught at draw time.
  expect_error(ap_pub_levels(factor(c("cancer", "normal")),
                             ap_pub_options(level_labels = list(cancer = "Normal"))),
               "same name")
})

test_that("the DA figures and legends use the level labels", {
  skip_pub()
  pub <- ap_pub_options(level_labels = list(b = "Tumour", a = "Control"))
  cc <- ap_da_concordance(fake_da())
  p <- ap_plot_da_primary(cc, "ancombc2", pub = pub)
  expect_equal(attr(p, "ap_heads")$text, c("Enriched in Tumour", "Depleted in Tumour"))
  expect_match(p$labels$x, "Tumour vs Control$")
  for (type in c("primary", "concordance", "volcano")) {
    leg <- ap_da_legend(cc, type, pub = pub)
    expect_equal(leg[1], "Contrast: Tumour vs Control (reference Control). 4 features tested in 20 samples.")
  }
  # Without labels the metadata values are used as they are.
  expect_equal(ap_da_legend(cc, "primary")[1],
               "Contrast: b vs a (reference a). 4 features tested in 20 samples.")
})

test_that("pairwise PERMANOVA legend lines use the level labels", {
  pw <- data.frame(metric = "bray_curtis", term = "dx", group1 = "adenoma", group2 = "cancer",
                   R2 = 0.01, pseudo_F = 2, p = 0.01, p_adj = 0.02, stringsAsFactors = FALSE)
  perm <- list(pairwise = pw)
  out <- ap_pairwise_lines(perm, "bray_curtis", "dx",
                           ap_pub_options(level_labels = list(cancer = "CRC")))
  expect_match(out[2], "^  adenoma vs CRC: ")
})


# Screen heatmaps ------------------------------------------------------------------------------

test_that("the publication screen draws one family per figure, every variable in one order", {
  skip_pub()
  f <- ap_fixture_screen_inputs()
  s <- ap_screen(alpha = f$alpha, beta = f$beta, permutations = 99L, n_resample = 20L)
  pa <- ap_plot_screen(s, publication = TRUE, family = "alpha")
  pb <- ap_plot_screen(s, publication = TRUE, family = "beta")
  expect_setequal(unique(pa$data$family), "alpha")
  expect_setequal(unique(pb$data$family), "beta")
  # Every test of the family is a tile; none is cut off.
  expect_equal(nrow(pa$data), sum(s$results$family == "alpha"))
  # Both figures list the variables in the same order.
  expect_identical(levels(pa$data$row), levels(pb$data$row))
  # Publication names, not Hill orders.
  expect_setequal(levels(pa$data$col), c("Richness", "Shannon"))
  # A dot on exactly the tests with q < 0.05.
  pts <- ggplot2::layer_data(pa, 2L)
  expect_equal(nrow(pts), sum(s$results$family == "alpha" & s$results$q < 0.05, na.rm = TRUE))
  # Negative adjusted effects are drawn as zero, never below.
  expect_true(all(pa$data$fill >= 0))
  expect_silent(ggplot2::ggplot_build(pa))
  expect_silent(ggplot2::ggplot_build(pb))
  expect_error(ap_plot_screen(s, publication = TRUE, family = "gamma"), "should be one of")
})

test_that("the screen legend explains the scale, stability and names for its family", {
  f <- ap_fixture_screen_inputs()
  s <- ap_screen(alpha = f$alpha, beta = f$beta, permutations = 99L, n_resample = 20L)
  la <- ap_screen_legend(s, "alpha")
  lb <- ap_screen_legend(s, "beta")
  expect_match(la[1], sprintf("alpha diversity: %d tests", sum(s$results$family == "alpha")))
  expect_match(la[1], sprintf("across all %d tests", nrow(s$results)))
  expect_true(any(grepl("not shared with the beta diversity figure", la, fixed = TRUE)))
  expect_true(any(grepl("Hill numbers of order 1 and 2", la, fixed = TRUE)))
  expect_false(any(grepl("Hill numbers", lb, fixed = TRUE)))
  expect_true(any(grepl("PERMANOVA R2 adjusted", lb, fixed = TRUE)))
  expect_match(la[length(la)], "^Hypothesis-generating")
})

test_that("publication and report use their own metric names", {
  expect_equal(ap_metric_label_pub(c("q0", "q1", "q2")), c("Richness", "Shannon", "Inverse Simpson"))
  expect_equal(ap_metric_label(c("q1", "q2")), c("Hill q1\n(Shannon)", "Hill q2\n(inverse Simpson)"))
})
