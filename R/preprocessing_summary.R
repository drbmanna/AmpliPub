# What each analysis was given.
#
# Tools with one preprocessing stage (MicrobiomeAnalyst filters and normalizes once, then
# runs everything on that table) are easy to describe because there is only one answer.
# AmpliPub gives each analysis the input its method assumes, so the answer differs by
# analysis and has to be written down, or a reader cannot tell which table a number came
# from. Every value is read from the results object; nothing is restated from defaults.

#' Which table, filter and normalization each analysis used
#'
#' One row per analysis in a workflow run: the input table, the filtering applied on
#' top of the shared feature table, and the normalization or transformation. It exists
#' because AmpliPub does not preprocess once for everything: diversity uses a rarefied
#' table, differential abundance uses raw counts and each method's own normalization,
#' and the composition figures use proportions or CLR.
#'
#' @param res The results list written by the workflow's analysis stage.
#' @return A data frame with columns `analysis`, `input`, `filter`, `normalization`.
#' @export
ap_preprocessing_summary <- function(res) {
  ap_assert(is.list(res), "`res` must be the workflow results list.")
  cfg <- res$config %||% list()
  an <- cfg$analysis %||% list()
  ft <- cfg$filter %||% list()
  seed <- an$seed %||% NA
  rows <- list()
  add <- function(analysis, input, filter, normalization) {
    rows[[length(rows) + 1L]] <<- data.frame(analysis = analysis, input = input,
                                             filter = filter, normalization = normalization,
                                             stringsAsFactors = FALSE)
  }
  pct <- function(v) sprintf("%g%%", 100 * v)

  shared <- c(
    if (nzchar(ft$include %||% "")) sprintf("taxa restricted to %s", ft$include),
    if (nzchar(ft$exclude %||% ""))
      sprintf("%s sequences removed", ap_and(trimws(strsplit(ft$exclude, ",")[[1]]))),
    if (nzchar(ft$drop_where %||% "")) sprintf("samples dropped where %s", ft$drop_where),
    if (isTRUE((ft$min_samples_fraction %||% 0) > 0))
      sprintf("features in fewer than %s of samples removed", pct(ft$min_samples_fraction)),
    if (isTRUE((ft$min_sample_reads %||% 0) > 0))
      sprintf("samples below %s reads removed", ft$min_sample_reads)
  )
  add("Feature table (every analysis below starts here)", "Filtered feature table",
      if (length(shared)) paste(shared, collapse = "; ") else "none",
      "none (counts)")

  rare <- if (!is.null(res$depth) && !is.na(res$depth)) {
    sprintf("rarefied once to %s reads in R, seed %s", format(res$depth, big.mark = ","), seed)
  } else "not rarefied"

  if (!is.null(res$alpha)) {
    add("Alpha diversity and its tests", "Rarefied table", "samples below the depth dropped", rare)
  }
  if (!is.null(res$alpha_repeated) && isTRUE(res$alpha_repeated$rarefied)) {
    add("Alpha diversity, repeated rarefaction (comparison)", "Filtered feature table",
        "samples below the depth dropped",
        sprintf("rarefied %d times to %s reads, values averaged, seed %s",
                as.integer(res$alpha_repeated$n_iter),
                format(res$alpha_repeated$depth, big.mark = ","), res$alpha_repeated$seed))
  }
  if (!is.null(res$beta)) {
    m <- res$beta$metrics
    norm <- rare
    if ("aitchison" %in% m) {
      norm <- paste0(norm, sprintf("; Aitchison: zeros replaced by %g, then CLR",
                                   res$beta$pseudocount %||% NA))
    }
    add(sprintf("Beta diversity (%s), PERMANOVA, dispersion, ordination",
                paste(m, collapse = ", ")),
        "Rarefied table", "samples below the depth dropped", norm)
  }
  if (!is.null(res$screen) || !is.null(res$explains)) {
    add("Metadata screen and competing-terms model", "Alpha and beta results above",
        "none beyond those", "as for alpha and beta diversity")
  }

  add("Composition bar charts", "Filtered feature table, collapsed to each rank",
      "top 15 taxa per group, combined; the rest pooled as Other",
      "proportions within each sample (total-sum scaling)")
  add("Composition heatmaps (publication)", "Filtered feature table, collapsed to each rank",
      "top 15 taxa per group, combined",
      "CLR on all taxa at the rank, zeros replaced by 0.5; group means centred on a baseline")

  if (!is.null(res$da)) {
    d <- res$da
    lv <- an$da_levels %||% character(0)
    filt <- c(if (length(lv)) sprintf("samples restricted to %s", paste(unlist(lv), collapse = ", ")),
              sprintf("features present in fewer than %s of those samples removed, once, for all methods",
                      pct(d$prv_cut)))
    method_norm <- c(
      ancombc2 = "ANCOM-BC2: bias-corrected log abundance from non-zero counts; calls must pass its pseudocount sensitivity analysis",
      aldex2 = "ALDEx2: CLR over Monte Carlo Dirichlet instances, all features as denominator",
      linda = "LinDA: CLR with LinDA's adaptive zero handling (package defaults)",
      maaslin2 = "MaAsLin2: total-sum scaling, then log"
    )
    add(sprintf("Differential abundance (%d features, %d samples)",
                as.integer(d$n_features), as.integer(d$n_samples)),
        "Filtered feature table, raw counts (never rarefied)",
        paste(filt, collapse = "; "),
        paste(method_norm[intersect(d$methods, names(method_norm))], collapse = ". "))
  }

  if (!is.null(res$normalization)) {
    ns <- res$normalization
    add("Normalization sensitivity", "Filtered feature table",
        "as the feature table",
        sprintf("each of %s separately, then PERMANOVA (Bray-Curtis; Aitchison for CLR), %d permutations",
                ap_and(ap_norm_short(ns$methods)), as.integer(ns$permutations)))
  }

  do.call(rbind, rows)
}

#' @keywords internal
ap_norm_short <- function(m) {
  lab <- c(tss = "TSS", css = "CSS", clr = "CLR", rarefy = "rarefying")
  out <- unname(lab[m])
  ifelse(is.na(out), m, out)
}
