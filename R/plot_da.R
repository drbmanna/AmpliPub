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
#' @param publication Draw the publication version: two columns of panels, each
#'   titled with its method and effect scale, the smallest q-values labelled in
#'   italics, calls that failed ANCOM-BC2's sensitivity analysis as open circles.
#'   Settings and counts go to [ap_da_legend()].
#' @param pub [ap_pub_options()] for `publication = TRUE`.
#'
#' @return A ggplot object.
#' @export
ap_plot_volcano <- function(da, contrast = NULL, label_top = 8L,
                            label_by = c("effect", "significance"),
                            publication = FALSE, pub = ap_pub_options()) {
  ap_assert(inherits(da, "ap_da"), "`da` must come from `ap_da()`, not {class(da)[1]}.")
  if (publication) return(ap_plot_volcano_pub(da, contrast, label_top, pub))
  label_by <- match.arg(label_by)

  df <- da$results
  contrast <- contrast %||% df$contrast[1]
  df <- df[df$contrast == contrast & !is.na(df$p_adj), , drop = FALSE]
  ap_assert(nrow(df) > 0L, "No results for contrast `{contrast}`.")

  df$neglog_q <- -log10(pmax(df$p_adj, .Machine$double.xmin))
  fs <- if ("failed_sensitivity" %in% names(df)) df$failed_sensitivity else FALSE
  df$call <- ifelse(fs, "failed sensitivity",
                    ifelse(!df$significant, "not significant",
                           ifelse(df$effect > 0, "enriched", "depleted")))
  df$call <- factor(df$call, levels = c("depleted", "not significant", "failed sensitivity",
                                        "enriched"))
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
      values = c(depleted = "#0072B2", `not significant` = "grey75",
                 `failed sensitivity` = "grey45", enriched = "#D55E00"),
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
  r <- r[r$contrast == concordance$contrast & r$feature %in% features & !is.na(r$effect), ,
         drop = FALSE]
  ap_assert(nrow(r) > 0L, "None of these features has an estimated effect to plot.")

  order_by <- stats::aggregate(effect ~ feature, data = r, FUN = function(v) mean(abs(v)))
  order_by <- order_by[order(-order_by$effect), ]
  keep <- utils::head(order_by$feature, max_features)
  r <- r[r$feature %in% keep, , drop = FALSE]

  ul <- ap_da_feature_labels(keep, r$taxon_label[match(keep, r$feature)])
  r$label <- ul[match(r$feature, keep)]
  label_order <- ul
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
#' @param publication Draw the publication version: one dot per method that
#'   called a feature, coloured by direction and sized by -log10 q, a small grey
#'   dot where the method tested it and found nothing, and an open circle where
#'   an ANCOM-BC2 call failed its sensitivity analysis. Features sharing a genus
#'   are told apart by a short ASV ID. Settings go to [ap_da_legend()].
#' @param pub [ap_pub_options()] for `publication = TRUE`.
#'
#' @return A ggplot object.
#' @export
ap_plot_concordance <- function(concordance, max_features = 40L, publication = FALSE,
                                pub = ap_pub_options()) {
  ap_assert(inherits(concordance, "ap_da_concordance"),
            "`concordance` must come from `ap_da_concordance()`, not {class(concordance)[1]}.")
  if (publication) return(ap_plot_concordance_pub(concordance, max_features, pub))
  feat <- concordance$features[concordance$features$n_methods > 0L, , drop = FALSE]
  ap_assert(nrow(feat) > 0L, "No feature was called by any method.")
  feat <- utils::head(feat, max_features)

  methods <- concordance$methods
  long <- do.call(rbind, lapply(methods, function(m) {
    data.frame(feature = feat$feature,
               label = ap_da_feature_labels(feat$feature, feat$taxon_label),
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

# Feature labels ------------------------------------------------------------------------------

# Differential abundance runs on ASVs, and several ASVs can share a genus (Baxter has three
# Prevotella ASVs). Keyed by genus alone, their rows land on top of each other. A label is
# made unique by adding the first characters of the feature ID, only where two shown
# features would otherwise share one.
#' @keywords internal
ap_da_feature_labels <- function(feature, taxon_label = NULL, id_chars = 6L) {
  feature <- as.character(feature)
  lab <- if (is.null(taxon_label)) feature else ifelse(is.na(taxon_label), feature, taxon_label)
  dup <- lab %in% lab[duplicated(lab)]
  lab[dup] <- sprintf("%s (ASV %s)", lab[dup], substr(feature[dup], 1L, id_chars))
  ap_assert(!anyDuplicated(lab),
            "Feature labels are still not unique after adding {id_chars} ID characters. Increase `id_chars`.")
  lab
}

# The same labels as plotmath: rank prefix removed, the name in italics, and the
# "(unassigned genus)" and "(ASV ...)" notes upright.
#' @keywords internal
ap_pub_feature_expr <- function(labels) {
  shown <- sub("^[a-z]__", "", labels)
  out <- lapply(shown, function(s) {
    notes <- regmatches(s, gregexpr(" \\([^()]*\\)", s))[[1]]
    name <- sub(" \\(.*$", "", s)
    if (length(notes) == 0L) return(bquote(italic(.(name))))
    note <- trimws(paste(notes, collapse = ""))
    bquote(italic(.(name)) ~ .(note))
  })
  stats::setNames(as.expression(out), labels)
}

#' @keywords internal
ap_da_method_title <- function(m) {
  lab <- c(ancombc2 = "ANCOM-BC2", aldex2 = "ALDEx2", linda = "LinDA", maaslin2 = "MaAsLin2")
  out <- unname(lab[m])
  ifelse(is.na(out), m, out)
}

# Publication concordance ----------------------------------------------------------------------

#' @keywords internal
ap_plot_concordance_pub <- function(concordance, max_features, pub) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  da <- concordance$da
  ap_assert(inherits(da, "ap_da"), "This concordance object does not carry its `ap_da` result.")
  r <- da$results[da$results$contrast == concordance$contrast, , drop = FALSE]
  if (!"failed_sensitivity" %in% names(r)) r$failed_sensitivity <- FALSE
  if (!"structural_zero" %in% names(r)) r$structural_zero <- FALSE
  if (!"direction" %in% names(r)) r$direction <- sign(r$effect)

  # Rows are features at least one method called. A failed-sensitivity ANCOM-BC2 result is
  # drawn only on those rows; as rows of its own it would bury the calls (145 on Baxter).
  shown <- unique(r$feature[r$significant])
  ap_assert(length(shown) > 0L, "No feature was called by any method, so there is nothing to draw.")
  feat <- concordance$features[concordance$features$feature %in% shown, , drop = FALSE]
  best_q <- suppressWarnings(tapply(r$p_adj[r$feature %in% shown], r$feature[r$feature %in% shown],
                                   min, na.rm = TRUE))
  best_q[!is.finite(best_q)] <- 0
  feat$best_q <- best_q[feat$feature]
  feat <- feat[order(-feat$n_methods, feat$best_q), , drop = FALSE]
  feat <- utils::head(feat, max_features)

  methods <- concordance$methods
  d <- r[r$feature %in% feat$feature & r$method %in% methods, , drop = FALSE]
  labs <- ap_da_feature_labels(feat$feature, feat$taxon_label)
  d$label <- factor(labs[match(d$feature, feat$feature)], levels = rev(labs))
  d$method <- factor(ap_da_method_title(d$method), levels = ap_da_method_title(methods))
  # Colour: direction. Shape: how the call was made (tested, absent from a group, or a
  # test that failed ANCOM-BC2's sensitivity analysis). Size: -log10 q, tested calls only.
  d$colour <- ifelse(d$failed_sensitivity, "Failed sensitivity analysis",
                     ifelse(d$direction > 0, "Enriched", "Depleted"))
  d$kind <- ifelse(d$structural_zero, "Absent from one group",
                   ifelse(d$failed_sensitivity, "Failed sensitivity analysis", "Tested"))
  d$neglog_q <- -log10(pmax(d$p_adj, .Machine$double.xmin))
  expr <- ap_pub_feature_expr(labs)

  tested <- d[(d$significant | d$failed_sensitivity) & !d$structural_zero, , drop = FALSE]
  absent <- d[d$structural_zero, , drop = FALSE]
  quiet <- d[!d$significant & !d$failed_sensitivity, , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$method, y = .data$label)) +
    ggplot2::geom_point(data = quiet, colour = "grey80", size = 0.6) +
    ggplot2::geom_point(data = tested,
                        ggplot2::aes(size = .data$neglog_q, colour = .data$colour, shape = .data$kind),
                        stroke = 0.5) +
    # Only a layer that draws something may draw a legend key; two layers keying the same
    # shape drew concentric circles.
    ggplot2::geom_point(data = absent,
                        ggplot2::aes(colour = .data$colour, shape = .data$kind),
                        size = 2.2, stroke = 0.5, show.legend = nrow(absent) > 0L) +
    ggplot2::scale_colour_manual(
      values = c(Enriched = "#B2182B", Depleted = "#2166AC",
                 `Failed sensitivity analysis` = "grey45"),
      breaks = c("Enriched", "Depleted"), name = NULL) +
    # A filled dot is an ordinary call and needs no key; only the exceptions are keyed.
    ggplot2::scale_shape_manual(values = c(Tested = 16, `Absent from one group` = 17,
                                           `Failed sensitivity analysis` = 1),
                                breaks = intersect(c("Absent from one group",
                                                     "Failed sensitivity analysis"), d$kind),
                                drop = TRUE, name = NULL) +
    ggplot2::guides(shape = ggplot2::guide_legend(
      override.aes = list(size = 2, colour = "grey45"))) +
    ggplot2::scale_size_continuous(name = expression(-log[10] ~ italic(q)), range = c(1, 3.5)) +
    ggplot2::scale_y_discrete(labels = expr[levels(d$label)]) +
    ggplot2::scale_x_discrete(drop = FALSE, position = "top") +
    ggplot2::labs(x = NULL, y = NULL) +
    ap_theme_pub(legend = "bottom") +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.2),
                   legend.box = "vertical", legend.spacing.y = ggplot2::unit(0, "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.text = ggplot2::element_text(size = 6),
                   # At 7 pt "ANCOM-BC2" and "ALDEx2" run into each other in 89 mm.
                   axis.text.x.top = ggplot2::element_text(size = 6))
  attr(p, "ap_pub_size") <- list(width = "single", height = min(247, max(60, 4.2 * nrow(feat) + 38)))
  p
}

# Publication volcano ----------------------------------------------------------------------------

#' @keywords internal
ap_plot_volcano_pub <- function(da, contrast, label_top, pub) {
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  df <- da$results
  contrast <- contrast %||% df$contrast[1]
  df <- df[df$contrast == contrast & !is.na(df$p_adj), , drop = FALSE]
  ap_assert(nrow(df) > 0L, "No results for contrast `{contrast}`.")
  if (!"failed_sensitivity" %in% names(df)) df$failed_sensitivity <- FALSE
  df$neglog_q <- -log10(pmax(df$p_adj, .Machine$double.xmin))
  df$call <- ifelse(df$failed_sensitivity, "Failed sensitivity analysis",
                    ifelse(!df$significant, "Not significant",
                           ifelse(df$effect > 0, "Enriched", "Depleted")))
  df$call <- factor(df$call, levels = c("Enriched", "Depleted", "Not significant",
                                        "Failed sensitivity analysis"))
  # "log fold change (natural log)" would nest brackets inside the panel title's own.
  scales <- gsub(" [(]([^)]*)[)]", ", \\1", tapply(df$effect_scale, df$method, `[`, 1))
  df$panel <- factor(sprintf("%s\n(%s)", ap_da_method_title(df$method), scales[df$method]),
                     levels = sprintf("%s\n(%s)", ap_da_method_title(da$methods), scales[da$methods]))

  sig <- df[df$significant, , drop = FALSE]
  # Labels are made unique per feature, not per row: a feature called by three
  # methods is one feature, not three sharing a name.
  uf <- unique(sig$feature)
  ul <- ap_da_feature_labels(uf, sig$taxon_label[match(uf, sig$feature)])
  sig$label <- ul[match(sig$feature, uf)]
  lab <- do.call(rbind, lapply(split(sig, sig$method), function(d) {
    utils::head(d[order(d$p_adj), , drop = FALSE], label_top)
  }))
  if (!is.null(lab) && nrow(lab) > 0L) {
    ex <- ap_pub_feature_expr(lab$label)
    lab$plotmath <- vapply(ex, function(e) paste(deparse(e), collapse = ""), character(1))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$effect, y = .data$neglog_q)) +
    ggplot2::geom_hline(yintercept = -log10(da$alpha), linetype = "dashed",
                        linewidth = 0.25, colour = "grey40") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$call, shape = .data$call),
                        size = 0.7, stroke = 0.4) +
    ggplot2::scale_colour_manual(values = c(Enriched = "#B2182B", Depleted = "#2166AC",
                                            `Not significant` = "grey75",
                                            `Failed sensitivity analysis` = "grey35"),
                                 drop = TRUE, name = NULL) +
    ggplot2::scale_shape_manual(values = c(Enriched = 16, Depleted = 16,
                                           `Not significant` = 16,
                                           `Failed sensitivity analysis` = 1),
                                drop = TRUE, name = NULL) +
    ggplot2::facet_wrap(~ panel, ncol = 2, scales = "free_x") +
    ggplot2::labs(x = "Effect size, on each method's own scale",
                  y = expression(-log[10] ~ "adjusted" ~ italic(p))) +
    ap_theme_pub(legend = "bottom") +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7, lineheight = 0.9),
                   legend.text = ggplot2::element_text(size = 6))
  if (!is.null(lab) && nrow(lab) > 0L && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = .data$plotmath),
                                      parse = TRUE, size = 6 / ggplot2::.pt, colour = "black",
                                      min.segment.length = 0, segment.size = 0.2,
                                      max.overlaps = Inf, seed = 1)
  }
  attr(p, "ap_pub_size") <- list(width = "double", height = 120)
  p
}

# Primary-method figures ----------------------------------------------------------------------

#' Main-text differential abundance figure from one primary method
#'
#' The features one method called, as that method's effect with its 95% interval
#' (`type = "effect"`, the main-text figure): enriched features above depleted, each block
#' sorted by effect. Agreement between methods is not repeated here; it is the concordance
#' figure's job ([ap_plot_concordance()]).
#'
#' `type = "abundance"` draws the per-sample relative abundance behind a few chosen calls
#' (`features`), in the same row order. It is for the handful of taxa a paper discusses in
#' the text: with dozens of rows it cannot be read, and for rare taxa most samples are zero,
#' so the detection count printed beside each group carries more than the box.
#'
#' The primary method is meant to be chosen before the results are seen (the workflow's
#' `da_primary`). Chosen afterwards, from whichever method called the most, it is the
#' selection the other methods are there to guard against.
#'
#' ALDEx2 cannot be the primary method for `type = "effect"`: the column AmpliPub stores as
#' its standard error is `diff.win`, a within-condition dispersion, and an interval built from
#' it would not be a confidence interval.
#'
#' @param concordance An `ap_da_concordance` result from [ap_da_concordance()].
#' @param primary The method whose calls are drawn. Default `"ancombc2"`.
#' @param type `"effect"` (default) or `"abundance"`.
#' @param x For `type = "abundance"`: the TreeSummarizedExperiment [ap_da()] ran on.
#' @param max_features Cap on rows in each direction, largest effects first. Default
#'   `25`: at most 25 enriched and 25 depleted features are drawn. Every call is in the
#'   results table, and the legend says how many are not shown.
#' @param features Feature IDs to draw, from the primary method's calls. Required for
#'   `type = "abundance"`, at most 8; optional for `type = "effect"`.
#' @param pub [ap_pub_options()].
#'
#' @return A ggplot object.
#' @export
ap_plot_da_primary <- function(concordance, primary = "ancombc2",
                               type = c("effect", "abundance"), x = NULL,
                               max_features = 25L, features = NULL, pub = ap_pub_options()) {
  type <- match.arg(type)
  ap_assert(inherits(pub, "ap_pub_options"), "`pub` must come from `ap_pub_options()`.")
  if (type == "abundance") {
    ap_assert(length(features) >= 1L && length(features) <= 8L, paste0(
      "`type = \"abundance\"` draws a few chosen calls: pass 1 to 8 feature IDs in `features`. ",
      "For every call, use the effect figure."))
  }
  rows <- ap_da_primary_rows(concordance, primary, max_features, features, pub)
  if (type == "effect") ap_plot_da_primary_effect(rows, concordance, primary, pub)
  else ap_plot_da_primary_abundance(rows, concordance, x, pub)
}

# One row per feature the primary method called, in drawing order (largest enriched effect
# first, largest depleted last). Structural zeros have no effect, so they sort to the
# extreme of their direction.
#' @keywords internal
ap_da_primary_rows <- function(concordance, primary, max_features, features = NULL, pub = NULL) {
  ap_assert(inherits(concordance, "ap_da_concordance"),
            "`concordance` must come from `ap_da_concordance()`, not {class(concordance)[1]}.")
  da <- concordance$da
  ap_assert(inherits(da, "ap_da"), "This concordance object does not carry its `ap_da` result.")
  ap_assert(length(primary) == 1L && primary %in% concordance$methods,
            "`primary` must be one of the methods that ran: {paste(concordance$methods, collapse = ', ')}.")
  r <- da$results[da$results$contrast == concordance$contrast, , drop = FALSE]
  if (!"structural_zero" %in% names(r)) r$structural_zero <- FALSE
  if (!"direction" %in% names(r)) r$direction <- sign(r$effect)
  p <- r[r$method == primary & r$significant, , drop = FALSE]
  ap_assert(nrow(p) > 0L, paste0(
    ap_da_method_title(primary), " called no feature at ", da$p_adj_method, "-adjusted p < ",
    da$alpha, ", so there is no primary-method figure to draw. Report that result, with the ",
    "concordance figure for what the other methods found."))
  if (!is.null(features)) {
    bad <- setdiff(features, p$feature)
    ap_assert(length(bad) == 0L, paste0(
      "{length(bad)} of `features` {?was/were} not called by ", ap_da_method_title(primary),
      ": {paste(bad, collapse = ', ')}."))
    p <- p[p$feature %in% features, , drop = FALSE]
  }

  key <-ifelse(p$structural_zero, p$direction * Inf, p$effect)
  p <- p[order(-key), , drop = FALSE]
  # At most `max_features` per direction, the largest effects (structural zeros first).
  up <- p$direction > 0
  n_up <- sum(up)
  n_down <- sum(!up)
  size <- ifelse(p$structural_zero, Inf, abs(p$effect))
  keep <- c(utils::head(which(up)[order(-size[up])], max_features),
            utils::head(which(!up)[order(-size[!up])], max_features))
  p <- p[sort(keep), , drop = FALSE]
  p$n_called_up <- n_up
  p$n_called_down <- n_down
  p$label <- ap_da_feature_labels(p$feature, p$taxon_label)
  # Enriched rows above depleted, each block topped by an empty row for its heading.
  up <- p$direction > 0
  nd <- sum(!up)
  p$y <- NA_real_
  p$y[!up] <- rev(seq_len(nd))
  p$y[up] <- rev(seq_len(sum(up))) + nd + if (nd > 0L) 1 else 0
  lev <- ap_pub_level_text(sub("_vs_.*$", "", concordance$contrast), pub)
  # Each heading sits on the side its points are on: enriched right, depleted left.
  heads <- data.frame(
    y = c(if (any(up)) max(p$y[up]) + 1, if (nd > 0L) nd + 1),
    text = c(if (any(up)) paste("Enriched in", lev), if (nd > 0L) paste("Depleted in", lev)),
    x = c(if (any(up)) Inf, if (nd > 0L) -Inf),
    hjust = c(if (any(up)) 1, if (nd > 0L) 0)
  )
  attr(p, "heads") <- heads
  p
}

#' @keywords internal
ap_plot_da_primary_effect <- function(rows, concordance, primary, pub = NULL) {
  ap_assert(primary != "aldex2", paste0(
    "ALDEx2 cannot be the primary method of an effect figure: its stored `se` is diff.win, a ",
    "within-condition dispersion, not a standard error. Use type = \"abundance\", or another method."))
  da <- concordance$da
  z <- stats::qnorm(0.975)
  rows$lo <- rows$effect - z * rows$se
  rows$hi <- rows$effect + z * rows$se
  rows$call <- ifelse(rows$direction > 0, "Enriched", "Depleted")
  tested <- rows[!rows$structural_zero, , drop = FALSE]
  ap_assert(all(is.finite(tested$se)),
            "{ap_da_method_title(primary)} returned no standard error for a called feature, so no interval can be drawn.")
  edge <- max(abs(c(tested$lo, tested$hi, 0)), na.rm = TRUE) * 1.08
  absent <- rows[rows$structural_zero, , drop = FALSE]
  if (nrow(absent)) absent$effect <- absent$direction * edge

  expr <- ap_pub_feature_expr(rows$label)
  lev <- ap_pub_level_text(sub("_vs_.*$", "", concordance$contrast), pub)
  scale <- tested$effect_scale[1] %||% rows$effect_scale[1]
  xlab <- paste0(toupper(substr(scale, 1, 1)), substring(scale, 2), ", ", lev, " vs ",
                 ap_pub_level_text(da$reference, pub))

  heads <- attr(rows, "heads")

  p <- ggplot2::ggplot(rows, ggplot2::aes(y = .data$y)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    ggplot2::geom_text(data = heads, ggplot2::aes(x = .data$x, y = .data$y, label = .data$text,
                                                  hjust = .data$hjust),
                       size = 6 / ggplot2::.pt, fontface = "bold", colour = "black") +
    ggplot2::geom_errorbar(data = tested, ggplot2::aes(xmin = .data$lo, xmax = .data$hi,
                                                       colour = .data$call),
                           orientation = "y", width = 0, linewidth = 0.45) +
    ggplot2::geom_point(data = tested, ggplot2::aes(x = .data$effect, colour = .data$call),
                        size = 1.6, shape = 16) +
    ggplot2::geom_point(data = absent, ggplot2::aes(x = .data$effect, colour = .data$call),
                        size = 1.8, shape = 17) +
    ggplot2::scale_colour_manual(values = c(Enriched = "#B2182B", Depleted = "#2166AC"),
                                 guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = rows$y, labels = expr, expand = ggplot2::expansion(add = 0.6),
      sec.axis = ggplot2::dup_axis(labels = ap_pub_q_expr(rows$p_adj),
                                   name = expression(italic(q)))) +
    ggplot2::labs(x = xlab, y = NULL) +
    ap_theme_pub() +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 6),
                   axis.text.y.right = ggplot2::element_text(size = 6, hjust = 0),
                   axis.title.y.right = ggplot2::element_text(size = 7, angle = 0, vjust = 1),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.2))
  # Long labels (unassigned genus, ASV suffix) take the width the data needs; go to 1.5
  # columns before the panel gets too narrow to read.
  width <- if (max(nchar(sub("^[a-z]__", "", rows$label))) > 26L) "onehalf" else "single"
  n_lines <- nrow(rows) + nrow(heads)
  attr(p, "ap_pub_size") <- list(width = width, height = min(247, max(22, 3 * n_lines + 18)))
  attr(p, "ap_n_called") <- c(up = rows$n_called_up[1], down = rows$n_called_down[1])
  attr(p, "ap_heads") <- heads
  p
}

# Adjusted p-values as plotmath, two significant figures: 0.012 as written, below 0.001 as
# 1.2 x 10^-6. A structural zero has no test and gets an empty label.
#' @keywords internal
ap_pub_q_expr <- function(q) {
  out <- lapply(q, function(v) {
    if (is.na(v)) return("")
    if (v >= 0.001) return(formatC(signif(v, 2), format = "fg", digits = 2))
    e <- floor(log10(v))
    m <- formatC(signif(v / 10^e, 2), format = "f", digits = 1)
    if (m == "10.0") { m <- "1.0"; e <- e + 1 }
    bquote(.(m) %*% 10^.(as.integer(e)))
  })
  as.expression(out)
}

#' @keywords internal
ap_plot_da_primary_abundance <- function(rows, concordance, x, pub) {
  ap_assert(!is.null(x), "`type = \"abundance\"` needs `x`, the object `ap_da()` ran on.")
  ap_validate(x)
  da <- concordance$da
  counts <- as.matrix(SummarizedExperiment::assay(x, "counts"))
  missing <- setdiff(rows$feature, rownames(counts))
  ap_assert(length(missing) == 0L,
            "{length(missing)} called feature{?s} {?is/are} not in `x`. Pass the object `ap_da()` ran on.")
  g <- SummarizedExperiment::colData(x)[[da$group]]
  keep <- !is.na(g) & g %in% c(da$reference, sub("_vs_.*$", "", concordance$contrast))
  ap_assert(sum(keep) == da$n_samples, paste0(
    "`x` has ", sum(keep), " samples in this contrast, but the differential abundance ran on ",
    da$n_samples, ". Pass the object `ap_da()` ran on."))
  counts <- counts[, keep, drop = FALSE]
  g <- droplevels(factor(g[keep], levels = c(da$reference, setdiff(unique(g[keep]), da$reference))))
  depth <- colSums(counts)
  rel <- 100 * sweep(counts[rows$feature, , drop = FALSE], 2, depth, "/")

  lv <- ap_pub_levels(g, pub)
  k <- nlevels(lv)
  off <- if (k == 1L) 0 else seq(-0.2, 0.2, length.out = k)
  long <- data.frame(
    feature = rep(rows$feature, times = ncol(rel)),
    value = as.vector(rel),
    grp = rep(lv, each = nrow(rel)),
    stringsAsFactors = FALSE
  )
  long$y <- rows$y[match(long$feature, rows$feature)] + off[as.integer(long$grp)]
  long$box <- interaction(long$feature, long$grp, drop = TRUE)

  # Below one read at the median depth the axis is linear, above it logarithmic, so the
  # zeros stay on the plot instead of being dropped by a log axis.
  sigma <- 100 / stats::median(depth)
  # Breaks inside the linear stretch would print on top of the 0.
  brks <- c(0, 10^(-3:2))
  brks <- brks[brks == 0 | (brks >= 3 * sigma & brks <= max(long$value) * 1.5)]
  expr <- ap_pub_feature_expr(rows$label)

  # Most samples of a rare taxon are zero and pile onto one point, so each box carries how
  # many samples the feature was detected in.
  det <- stats::aggregate(value ~ feature + grp, data = long,
                          FUN = function(v) sprintf("%d/%d", sum(v > 0), length(v)))
  det$y <- rows$y[match(det$feature, rows$feature)] + off[as.integer(det$grp)]

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$y)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = .data$box, fill = .data$grp), orientation = "y",
                          width = 0.8 / k * 0.85, linewidth = 0.3, colour = "black",
                          outlier.shape = NA) +
    ggplot2::geom_point(size = 0.6, alpha = 0.5, colour = "black", stroke = 0,
                        position = ggplot2::position_jitter(width = 0, height = 0.06, seed = 1)) +
    ggplot2::geom_text(data = det, ggplot2::aes(x = Inf, y = .data$y, label = .data$value),
                       hjust = -0.15, size = 6 / ggplot2::.pt, colour = "black") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_fill_manual(values = ap_pub_palette(k, pub$palette), name = NULL,
                               breaks = rev(levels(lv))) +
    ggplot2::scale_x_continuous(trans = scales::pseudo_log_trans(sigma = sigma, base = 10),
                                breaks = brks, labels = function(v) format(v, drop0trailing = TRUE,
                                                                           scientific = FALSE)) +
    ggplot2::scale_y_continuous(breaks = rows$y, labels = expr,
                                expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::labs(x = "Relative abundance (%)", y = NULL) +
    ap_theme_pub(legend = "bottom") +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.2),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   plot.margin = ggplot2::margin(2, 10, 2, 2, unit = "mm"))
  attr(p, "ap_pub_size") <- list(width = "single", height = min(247, max(35, 7 * nrow(rows) + 26)))
  attr(p, "ap_sigma") <- sigma
  p
}

#' Legend text for a publication differential abundance figure
#'
#' The settings and counts a reader needs next to the concordance or volcano figure:
#' the contrast and reference, how many features and samples were tested, the filter and
#' correction, the count each method called, how many ANCOM-BC2 calls failed its
#' sensitivity analysis, the consensus, and the full IDs of features shown by short ID.
#'
#' @param concordance An `ap_da_concordance` result.
#' @param type `"concordance"`, `"volcano"`, or `"primary"` (the figure from
#'   [ap_plot_da_primary()]).
#' @param primary,max_features For `type = "primary"`: the same values given to
#'   [ap_plot_da_primary()].
#' @param pub [ap_pub_options()], for the group names (`level_labels`).
#' @return A character vector, one line per entry.
#' @export
ap_da_legend <- function(concordance, type = c("concordance", "volcano", "primary"),
                         primary = "ancombc2", max_features = 25L, pub = ap_pub_options()) {
  type <- match.arg(type)
  ap_assert(inherits(concordance, "ap_da_concordance"),
            "`concordance` must come from `ap_da_concordance()`.")
  if (type == "primary") return(ap_da_legend_primary(concordance, primary, max_features, pub))
  da <- concordance$da
  r <- da$results[da$results$contrast == concordance$contrast, , drop = FALSE]
  if (!"failed_sensitivity" %in% names(r)) r$failed_sensitivity <- FALSE
  n_by <- vapply(concordance$methods, function(m) sum(r$significant[r$method == m]), numeric(1))
  n_fail <- sum(r$failed_sensitivity)
  n_sz <- if ("structural_zero" %in% names(r)) sum(r$structural_zero) else 0L
  sig_ids <- unique(r$feature[r$significant])
  labs <- if (length(sig_ids)) {
    ap_da_feature_labels(sig_ids, r$taxon_label[match(sig_ids, r$feature)])
  } else character(0)
  short <- grepl("\\(ASV ", labs)
  c(
    ap_da_contrast_line(concordance, pub),
    sprintf(paste0("Features present in fewer than %g%% of samples were removed once, for all ",
                   "methods. %s-adjusted p < %g."), 100 * da$prv_cut, da$p_adj_method, da$alpha),
    paste0("Called: ", paste(sprintf("%s %d", ap_da_method_title(names(n_by)), as.integer(n_by)),
                             collapse = ", "), "."),
    if (n_fail > 0L) sprintf(paste0("%d ANCOM-BC2 call%s had adjusted p < %g but failed its ",
                                    "pseudocount sensitivity analysis and %s not counted%s."),
                             n_fail, if (n_fail == 1L) "" else "s", da$alpha,
                             if (n_fail == 1L) "is" else "are",
                             if (type == "concordance") " (open circles where the feature was called by another method)" else " (open circles)"),
    if (n_sz > 0L) sprintf(paste0("%d ANCOM-BC2 call%s %s structural zero%s: absent from one ",
                                  "group, declared differentially abundant without an effect ",
                                  "or p-value%s."),
                           n_sz, if (n_sz == 1L) "" else "s", if (n_sz == 1L) "is a" else "are",
                           if (n_sz == 1L) "" else "s",
                           if (type == "concordance") " (triangles)" else "; not placed on the volcano"),
    sprintf("Consensus (called by at least %d methods, same direction): %d feature%s.",
            as.integer(concordance$min_methods), length(concordance$consensus),
            if (length(concordance$consensus) == 1L) "" else "s"),
    if (type == "concordance") {
      paste0("Dot colour gives the direction relative to ", ap_pub_level_text(da$reference, pub),
             "; dot size gives -log10 q. Effect magnitudes are not compared across methods, ",
             "because each method reports a different quantity.")
    } else {
      paste0("Each panel uses its method's own effect scale; panels cannot be compared by ",
             "position. Dashed line: adjusted p = ", da$alpha, ".")
    },
    if (any(short)) paste0("Full feature IDs: ",
                           paste(sprintf("%s = %s", labs[short], sig_ids[short]), collapse = "; "),
                           ".")
  )
}

#' @keywords internal
ap_da_legend_primary <- function(concordance, primary, max_features, pub) {
  rows <- ap_da_primary_rows(concordance, primary, max_features, pub = pub)
  da <- concordance$da
  r <- da$results[da$results$contrast == concordance$contrast & da$results$method == primary, ,
                  drop = FALSE]
  n_fail <- if ("failed_sensitivity" %in% names(r)) sum(r$failed_sensitivity) else 0L
  n_up <- rows$n_called_up[1]
  n_down <- rows$n_called_down[1]
  s_up <- sum(rows$direction > 0)
  s_down <- sum(rows$direction < 0)
  lev <- ap_pub_level_text(sub("_vs_.*$", "", concordance$contrast), pub)
  scale <- r$effect_scale[!is.na(r$effect_scale)][1]
  short <- grepl("[(]ASV ", rows$label)
  n_sz <- sum(rows$structural_zero)
  c(
    ap_da_contrast_line(concordance, pub),
    sprintf(paste0("Features present in fewer than %g%% of samples were removed once, for all ",
                   "methods. %s-adjusted p < %g."), 100 * da$prv_cut, da$p_adj_method, da$alpha),
    sprintf(paste0("Primary method %s. Points: %s with 95%% interval (estimate %s 1.96 SE, not ",
                   "adjusted for multiple testing). Right: q, the %s-adjusted p-value from %s's ",
                   "own test."),
            ap_da_method_title(primary), scale, "\u00b1", da$p_adj_method,
            ap_da_method_title(primary)),
    if (primary == "ancombc2") {
      paste0("A feature counts as called only if it also passed ANCOM-BC2's pseudocount ",
             "sensitivity analysis",
             if (n_fail == 0L) "." else sprintf("; %d with adjusted p < %g did not and %s not shown.",
                                               as.integer(n_fail), da$alpha,
                                               if (n_fail == 1L) "is" else "are"))
    },
    if (s_up < n_up || s_down < n_down) {
      sprintf(paste0("Shown: the %d largest of %d enriched and %d largest of %d depleted calls in ",
                     "%s; every call is in da_results.tsv."), s_up, n_up, s_down, n_down, lev)
    } else if (n_up + n_down == 1L) {
      sprintf("One feature called (%s in %s).", if (n_up == 1L) "enriched" else "depleted", lev)
    } else {
      sprintf("All %d calls shown (%d enriched, %d depleted in %s).", n_up + n_down, n_up, n_down, lev)
    },
    if (n_sz > 0L) {
      sprintf(paste0("Triangles at the axis edge: %d feature%s absent from one group ",
                     "(structural zero), called without an effect estimate."),
              n_sz, if (n_sz == 1L) "" else "s")
    },
    "Calls by the other methods are compared in the concordance figure.",
    if (any(short)) paste0("Full feature IDs: ",
                           paste(sprintf("%s = %s", sub("^[a-z]__", "", rows$label[short]),
                                         rows$feature[short]),
                                 collapse = "; "), ".")
  )
}

# "Contrast: cancer vs normal (reference normal). 391 features tested in 292 samples."
#' @keywords internal
ap_da_contrast_line <- function(concordance, pub) {
  da <- concordance$da
  lev <- ap_pub_level_text(sub("_vs_.*$", "", concordance$contrast), pub)
  ref <- ap_pub_level_text(da$reference, pub)
  sprintf("Contrast: %s vs %s (reference %s). %d features tested in %d samples.",
          lev, ref, ref, as.integer(da$n_features), as.integer(da$n_samples))
}
