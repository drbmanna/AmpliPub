# Differential abundance.
#
# The most contested methodology in the field. Benchmarks disagree with each
# other, methods disagree with each other on the same data, and a paper that
# runs one method and reports its significant list is reporting a property of
# that method as much as a property of the biology. So AmpliPub runs four and
# reports where they agree.
#
# Three design decisions follow from that.
#
# 1. Prevalence filtering happens ONCE, here, before any method runs, and every
#    method is then told not to filter again. Each package ships its own default
#    (ancombc2 prv_cut = 0.1, linda prev.filter = 0, Maaslin2 min_prevalence =
#    0.1), and letting each apply its own means the methods test different
#    feature sets and their FDR denominators are not comparable. Comparing
#    q-values across different denominators is not a comparison.
#
# 2. Effect sizes are NOT put on a common scale. ANCOM-BC2 reports a natural-log
#    fold change of bias-corrected abundance, LinDA a log2 fold change, MaAsLin2
#    a coefficient on log-transformed relative abundance, and ALDEx2 a
#    standardised median CLR difference. These are different quantities, and
#    rescaling them into one column would invent a comparability that does not
#    exist. Each row carries its own `effect_scale`.
#
# 3. Concordance is therefore computed on direction and significance, which are
#    comparable, rather than on effect magnitude, which is not.
#
# Nothing here rarefies. Rarefaction discards data and distorts the variance
# structure these models estimate; it is defensible for richness and wrong here.

#' @keywords internal
ap_da_methods <- function() c("ancombc2", "aldex2", "linda", "maaslin2")

#' Differential abundance
#'
#' Runs one or more differential abundance methods on a common, pre-filtered
#' feature set and returns a single result frame.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param group Metadata variable defining the comparison.
#' @param method Methods to run. Default: all four.
#' @param covariates Optional variables to adjust for. Not supported by ALDEx2,
#'   which is skipped with a note when covariates are given.
#' @param subject Optional variable identifying repeated measures. Supported by
#'   ANCOM-BC2, LinDA and MaAsLin2 as a random intercept.
#' @param reference Level of `group` to compare against. Defaults to the first
#'   in sorted order, and is recorded so the sign of an effect is unambiguous.
#' @param prv_cut Minimum fraction of samples in which a feature must be present
#'   to be tested. Applied once, to all methods. Default `0.1`.
#' @param alpha Significance threshold on the adjusted p-value. Default `0.05`.
#' @param p_adj_method Multiple testing correction. Default `"BH"`.
#' @param seed Random seed, recorded with the result.
#' @param mc_samples ALDEx2 Monte Carlo instances. Default `128`.
#'
#' @return An object of class `ap_da`: a list with `results` (one row per
#'   feature per method) and the settings used.
#' @export
ap_da <- function(x,
                  group,
                  method = ap_da_methods(),
                  covariates = NULL,
                  subject = NULL,
                  reference = NULL,
                  prv_cut = 0.1,
                  alpha = 0.05,
                  p_adj_method = "BH",
                  seed = 1L,
                  mc_samples = 128L) {

  bad <- setdiff(method, ap_da_methods())
  ap_assert(length(bad) == 0L,
            "{cli::qty(length(bad))}Unknown method{?s}: {paste(bad, collapse = ', ')}. Known: {paste(ap_da_methods(), collapse = ', ')}.")

  meta <- as.data.frame(SummarizedExperiment::colData(x))
  ap_assert(group %in% names(meta),
            "Variable `{group}` is not in the sample metadata. Available: {paste(names(meta), collapse = ', ')}.")
  for (v in c(covariates, subject)) {
    ap_assert(v %in% names(meta), "Variable `{v}` is not in the sample metadata.")
  }

  counts <- SummarizedExperiment::assay(x, "counts")
  ap_check_count_matrix(counts)
  ap_assert(all(counts == round(counts)),
            paste0("Differential abundance methods here model counts. This table holds ",
                   "fractional values, so it has already been normalised."))

  # One filter, before any method runs.
  keep_samples <- rownames(meta)[stats::complete.cases(meta[, c(group, covariates, subject), drop = FALSE])]
  keep_samples <- intersect(colnames(counts), keep_samples)
  if (length(keep_samples) < ncol(counts)) {
    cli::cli_inform("{ncol(counts) - length(keep_samples)} sample{?s} dropped for missing metadata.")
  }
  counts <- counts[, keep_samples, drop = FALSE]
  meta <- meta[keep_samples, , drop = FALSE]

  prevalence <- rowMeans(counts > 0)
  keep <- prevalence >= prv_cut
  ap_assert(sum(keep) > 1L,
            "Prevalence filter at {prv_cut} leaves {sum(keep)} feature{?s}. Lower `prv_cut`.")
  cli::cli_inform(paste0(
    "Prevalence filter at {prv_cut}: {sum(keep)} of {length(keep)} features tested. ",
    "The same set goes to every method, so the FDR denominators are comparable."
  ))
  counts <- counts[keep, , drop = FALSE]

  g <- meta[[group]]
  ap_assert(!is.numeric(g) || length(unique(g)) <= 10L,
            "`{group}` is continuous. This layer compares groups; use a categorical variable.")
  g <- factor(as.character(g))
  reference <- reference %||% levels(g)[1]
  ap_assert(reference %in% levels(g),
            "Reference `{reference}` is not a level of `{group}`. Levels: {paste(levels(g), collapse = ', ')}.")
  g <- stats::relevel(g, ref = reference)
  meta[[group]] <- g

  if (nlevels(g) > 2L) {
    cli::cli_inform(paste0(
      "`{group}` has {nlevels(g)} levels. Each level is compared against `{reference}`, ",
      "so a feature may be called for one contrast and not another."
    ))
  }

  rows <- list()
  skipped <- character(0)
  reasons <- character(0)
  for (m in method) {
    res <- tryCatch(
      ap_da_run_one(m, counts, meta, group, covariates, subject, reference,
                    alpha, p_adj_method, seed, mc_samples),
      error = function(e) {
        ap_warn("Method {m} failed and was skipped: {conditionMessage(e)}")
        reasons[[m]] <<- conditionMessage(e)
        NULL
      }
    )
    if (is.null(res)) {
      skipped <- c(skipped, m)
    } else {
      # A backend that returns feature IDs which are not the ones it was given
      # produces a result that looks complete and silently fails to join to the
      # other methods, which then makes the concordance meaningless. Maaslin2
      # did exactly this via make.names(). Checked here so any future backend
      # has to fail loudly rather than quietly.
      stray <- setdiff(unique(res$feature), rownames(counts))
      ap_assert(
        length(stray) == 0L,
        paste0("Method {m} returned {length(stray)} feature ID{?s} that were not in the ",
               "table it was given, for example {paste(utils::head(stray, 3), collapse = ', ')}. ",
               "Its results cannot be joined to the other methods.")
      )
      rows[[length(rows) + 1L]] <- res
    }
  }
  # Carrying the reasons into the abort matters. "Every method failed" on its
  # own sends the reader back to run the methods one at a time to find out why.
  ap_assert(
    length(rows) > 0L,
    paste0("Every method failed, so there is nothing to report.\n",
           "{paste(sprintf('  %s: %s', names(reasons), reasons), collapse = '\n')}")
  )

  results <- do.call(rbind, rows)
  # Each package names its contrast differently. Normalising here, rather than
  # in each adapter, is what lets the concordance layer line the methods up at
  # all; without it every method looks like it tested a different comparison.
  results$contrast <- paste0(results$contrast, "_vs_", reference)
  results$prevalence <- prevalence[match(results$feature, names(prevalence))]
  results$significant <- !is.na(results$p_adj) & results$p_adj < alpha
  rownames(results) <- NULL

  # Taxonomy carried along, so a result is readable without a second join.
  rd <- as.data.frame(SummarizedExperiment::rowData(x))
  if ("genus" %in% colnames(rd)) {
    results$taxon_label <- ap_taxon_label(rd, "genus")[match(results$feature, rownames(rd))]
  }

  structure(
    list(results = results,
         group = group, reference = reference, levels = levels(g),
         methods = setdiff(method, skipped), skipped = skipped,
         covariates = covariates, subject = subject,
         prv_cut = prv_cut, alpha = alpha, p_adj_method = p_adj_method,
         seed = seed, n_features = nrow(counts), n_samples = ncol(counts)),
    class = "ap_da"
  )
}

#' @keywords internal
ap_da_run_one <- function(m, counts, meta, group, covariates, subject, reference,
                          alpha, p_adj_method, seed, mc_samples) {
  set.seed(seed)
  switch(
    m,
    ancombc2 = ap_da_ancombc2(counts, meta, group, covariates, subject, alpha, p_adj_method, seed),
    aldex2 = ap_da_aldex2(counts, meta, group, covariates, subject, reference, mc_samples),
    linda = ap_da_linda(counts, meta, group, covariates, subject, alpha, p_adj_method),
    maaslin2 = ap_da_maaslin2(counts, meta, group, covariates, subject, reference, p_adj_method)
  )
}

# The shape every backend returns. Effect sizes stay on their native scale and
# say what that scale is.
#' @keywords internal
ap_da_row <- function(feature, method, effect, effect_scale, se, statistic,
                      p, p_adj, contrast, note = NA_character_) {
  data.frame(
    feature = as.character(feature),
    method = method,
    contrast = contrast,
    effect = as.numeric(effect),
    effect_scale = effect_scale,
    se = as.numeric(se),
    statistic = as.numeric(statistic),
    p = as.numeric(p),
    p_adj = as.numeric(p_adj),
    note = note,
    stringsAsFactors = FALSE
  )
}

#' @export
print.ap_da <- function(x, ...) {
  cli::cli_h1("Differential abundance: {x$group}")
  cli::cli_text("Reference level: {.strong {x$reference}} | ",
                "{x$n_features} features x {x$n_samples} samples | ",
                "prevalence filter {x$prv_cut} | {x$p_adj_method}, alpha {x$alpha}")
  if (!is.null(x$covariates)) cli::cli_text("Adjusted for: {paste(x$covariates, collapse = ', ')}")
  if (!is.null(x$subject)) cli::cli_text("Random intercept per: {x$subject}")
  if (length(x$skipped) > 0L) {
    cli::cli_alert_warning("Skipped: {paste(x$skipped, collapse = ', ')}")
  }

  tab <- stats::aggregate(significant ~ method + contrast, data = x$results, FUN = sum)
  names(tab)[3] <- "n_significant"
  tab$n_tested <- stats::aggregate(significant ~ method + contrast,
                                   data = x$results, FUN = length)$significant
  print(tab, row.names = FALSE)

  cli::cli_text("")
  cli::cli_alert_info(paste0(
    "Counts differing between methods is the normal outcome, not a bug. ",
    "Use `ap_da_concordance()` for the consensus set, and report the ",
    "disagreement rather than the most generous method."
  ))
  invisible(x)
}
