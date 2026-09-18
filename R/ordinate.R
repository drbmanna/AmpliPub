# Ordination.
#
# Both methods here carry a number that belongs on the figure and is routinely
# left off: PCoA axes carry the proportion of variation each explains, and NMDS
# carries stress. A PCoA whose first axis explains 9% looks exactly like one
# explaining 60% unless the axis says so, and readers reasonably assume the
# latter. Stress says whether the NMDS layout can be trusted at all.

#' Ordinate a distance matrix
#'
#' Principal coordinates analysis (PCoA) or non-metric multidimensional scaling
#' (NMDS).
#'
#' PCoA places samples so that Euclidean distance in the plot approximates the
#' input distance, and reports how much variation each axis carries. Non-
#' Euclidean distances such as Bray-Curtis can produce negative eigenvalues;
#' the proportion they account for is reported rather than dropped, because a
#' large negative fraction means the projection is distorting the data.
#'
#' NMDS preserves rank order instead and reports stress. The usual reading is
#' that stress below 0.1 is good, 0.1 to 0.2 is usable, and above 0.2 means the
#' plot should not be interpreted as a map.
#'
#' @param beta An `ap_beta` object from [ap_beta()], or a `stats::dist`.
#' @param metric Which metric to ordinate when `beta` holds several.
#' @param method `"pcoa"` (default) or `"nmds"`.
#' @param k Number of dimensions. Default `2`.
#' @param seed Random seed for NMDS, recorded with the result.
#' @param trymax NMDS restarts. Default `50`.
#'
#' @return An object of class `ap_ordination`.
#' @export
ap_ordinate <- function(beta, metric = NULL, method = c("pcoa", "nmds"),
                        k = 2L, seed = 1L, trymax = 50L) {
  method <- match.arg(method)

  if (inherits(beta, "ap_beta")) {
    metric <- metric %||% names(beta$distances)[1]
    ap_assert(metric %in% names(beta$distances),
              "Metric `{metric}` is not in this ap_beta object. Available: {paste(names(beta$distances), collapse = ', ')}.")
    d <- beta$distances[[metric]]
    meta <- beta$metadata
  } else {
    ap_assert(inherits(beta, "dist"),
              "`beta` must be an `ap_beta` object or a `stats::dist`, not {class(beta)[1]}.")
    d <- beta
    metric <- metric %||% "distance"
    meta <- NULL
  }

  set.seed(seed)
  if (method == "pcoa") {
    pc <- stats::cmdscale(d, k = k, eig = TRUE)
    eig <- pc$eig
    # Non-Euclidean distances give negative eigenvalues. Denominator uses the
    # positive eigenvalues only, which is the standard convention, and the
    # negative fraction is reported so the distortion is visible.
    pos <- eig[eig > 0]
    prop <- eig / sum(pos)
    neg_frac <- sum(abs(eig[eig < 0])) / sum(abs(eig))
    coords <- pc$points
    colnames(coords) <- paste0("PCo", seq_len(ncol(coords)))
    out <- list(coords = coords, eig = eig, prop_explained = prop,
                negative_eigenvalue_fraction = neg_frac, stress = NA_real_)
  } else {
    ord <- vegan::metaMDS(d, k = k, trymax = trymax, trace = 0, autotransform = FALSE)
    coords <- vegan::scores(ord, display = "sites")
    colnames(coords) <- paste0("NMDS", seq_len(ncol(coords)))
    out <- list(coords = coords, eig = NULL, prop_explained = NULL,
                negative_eigenvalue_fraction = NA_real_, stress = ord$stress,
                converged = ord$converged)
  }

  structure(
    c(out, list(method = method, metric = metric, k = k, seed = seed,
                metadata = meta, n_samples = nrow(out$coords))),
    class = "ap_ordination"
  )
}

#' @keywords internal
ap_axis_label <- function(ord, axis) {
  nm <- colnames(ord$coords)[axis]
  if (ord$method == "pcoa") {
    sprintf("%s (%.1f%%)", nm, 100 * ord$prop_explained[axis])
  } else {
    nm
  }
}

#' Diagnostics for one or more ordinations
#'
#' An ordination always produces a picture. Whether that picture can be read as
#' a map is a separate question, and it is answered by the stress (NMDS) or by
#' the negative eigenvalue mass and the variation on the first two axes (PCoA).
#' Before 2026-09-18 these numbers were computed inside [ap_ordinate()] and
#' discarded, so every report showed ordination panels with no indication of
#' whether the projection was faithful.
#'
#' @param ords A single `ap_ordination` or a list of them.
#'
#' @return A data frame with one row per ordination: `metric`, `method`,
#'   `n_samples`, `stress`, `converged`, `negative_eigenvalue_fraction`,
#'   `axes_1_2_explained`, `verdict` and `interpretation`.
#' @export
ap_ordination_diagnostics <- function(ords) {
  if (inherits(ords, "ap_ordination")) ords <- list(ords)
  do.call(rbind, lapply(ords, ap_ordination_interpret))
}

# Thresholds. NMDS stress follows the conventional rule of thumb, below 0.1 a
# good representation and above 0.2 one that should not be read as a map. The
# 5% negative eigenvalue mass and the 20% on the first two axes are choices
# recorded with the result, not standards: they mark the point where saying
# "the groups sit apart in this plot" stops being supportable.
#' @keywords internal
ap_ordination_interpret <- function(x) {
  base <- data.frame(
    metric = x$metric, method = x$method, n_samples = x$n_samples,
    stress = NA_real_, converged = NA,
    negative_eigenvalue_fraction = NA_real_, axes_1_2_explained = NA_real_,
    verdict = NA_character_, interpretation = NA_character_,
    stringsAsFactors = FALSE
  )

  if (identical(x$method, "pcoa")) {
    neg <- x$negative_eigenvalue_fraction
    ax12 <- sum(x$prop_explained[1:2])
    base$negative_eigenvalue_fraction <- neg
    base$axes_1_2_explained <- ax12
    distorted <- neg > 0.05
    flat <- ax12 < 0.20

    if (distorted && flat) {
      base$verdict <- "distorted and low-variance"
      base$interpretation <- sprintf(
        paste0("%.1f%% of the eigenvalue mass is negative, so this distance is not ",
               "Euclidean and the projection distorts it, and the first two axes carry ",
               "only %.1f%% of the variation. Distances read off this plot are not ",
               "reliable. Use the PERMANOVA result rather than the picture, and ",
               "consider NMDS, which makes no Euclidean assumption."),
        100 * neg, 100 * ax12)
    } else if (distorted) {
      base$verdict <- "distorted projection"
      base$interpretation <- sprintf(
        paste0("%.1f%% of the eigenvalue mass is negative. This distance is not ",
               "Euclidean, so the projection distorts it. The first two axes carry ",
               "%.1f%% of the variation. Consider NMDS, which makes no Euclidean ",
               "assumption."),
        100 * neg, 100 * ax12)
    } else if (flat) {
      base$verdict <- "low variance on the plotted axes"
      base$interpretation <- sprintf(
        paste0("The projection is close to Euclidean (%.1f%% negative eigenvalue mass), ",
               "but the first two axes carry only %.1f%% of the variation. Most of the ",
               "structure is in axes that are not plotted, so groups that look separated ",
               "here may not be, and groups that overlap here may still differ."),
        100 * neg, 100 * ax12)
    } else {
      base$verdict <- "readable"
      base$interpretation <- sprintf(
        paste0("The projection is close to Euclidean (%.1f%% negative eigenvalue mass) ",
               "and the first two axes carry %.1f%% of the variation. Distances in this ",
               "plot can be read as approximate distances in the full space."),
        100 * neg, 100 * ax12)
    }
    return(base)
  }

  stress <- x$stress
  converged <- isTRUE(x$converged)
  base$stress <- stress
  base$converged <- converged

  if (is.na(stress)) {
    base$verdict <- "not available"
    base$interpretation <- "NMDS returned no stress value, so the fit cannot be judged."
    return(base)
  }

  if (stress >= 0.2) {
    base$verdict <- "too high to read as a map"
    base$interpretation <- sprintf(
      paste0("Stress is %.3f. Above 0.2 the ordination is not a faithful summary of the ",
             "distances, so the arrangement of points should not be interpreted. Report ",
             "the PERMANOVA and the distances, not this plot."),
      stress)
  } else if (stress >= 0.1) {
    base$verdict <- "usable"
    base$interpretation <- sprintf(
      paste0("Stress is %.3f. Between 0.1 and 0.2 the map is usable for broad pattern ",
             "but not for fine distinctions. Read large separations, not small ones."),
      stress)
  } else {
    base$verdict <- "good"
    base$interpretation <- sprintf(
      paste0("Stress is %.3f. Below 0.1 the two-dimensional arrangement represents the ",
             "distances well, so relative positions can be read directly."),
      stress)
  }

  if (!converged) {
    base$verdict <- paste0(base$verdict, ", did not converge")
    base$interpretation <- paste0(
      base$interpretation,
      " NMDS did not converge within `trymax` attempts, so this solution may not be the ",
      "best available. Raise `trymax` and rerun before reporting it.")
  }
  base
}

#' @export
print.ap_ordination <- function(x, ...) {
  cli::cli_h1("Ordination: {toupper(x$method)} on {x$metric}")
  cli::cli_text("{x$n_samples} samples in {x$k} dimension{?s}, seed {x$seed}")
  if (x$method == "pcoa") {
    shown <- seq_len(min(5L, length(x$prop_explained)))
    cli::cli_text("Variation explained: {paste(sprintf('PCo%d %.1f%%', shown, 100 * x$prop_explained[shown]), collapse = ', ')}")
    cli::cli_text("First two axes together: {sprintf('%.1f%%', 100 * sum(x$prop_explained[1:2]))}")
    if (x$negative_eigenvalue_fraction > 0.05) {
      cli::cli_alert_warning(paste0(
        "{sprintf('%.1f%%', 100 * x$negative_eigenvalue_fraction)} of the eigenvalue mass ",
        "is negative. This distance is not Euclidean, so the projection distorts it. ",
        "Consider NMDS, which makes no Euclidean assumption."
      ))
    }
  } else {
    cli::cli_text("Stress: {sprintf('%.4f', x$stress)}")
    verdict <- if (x$stress < 0.1) "good" else if (x$stress < 0.2) "usable" else
      "too high to read as a map"
    cli::cli_text("Interpretation: {verdict}")
    if (isFALSE(x$converged)) cli::cli_alert_warning("NMDS did not converge; raise `trymax`.")
  }
  invisible(x)
}
