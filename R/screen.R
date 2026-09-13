# An exploratory screen across every metadata variable.
#
# "Run every test and show which variable explains the data best" is a fair
# request with a trap in it. Rank by p and the sample size picks the winner: at
# n = 487 a variable explaining 1% of the variation gets p = 0.001, and so does
# a batch proxy explaining 1.2%. Correct per test instead of across the screen
# and several hundred tests hand back a couple of dozen false positives.
#
# So the screen ranks by effect size on one scale, variance explained, applies
# one BH correction across everything it ran, prints how many tests that was,
# and reports how stable each rank is under subsampling. It is labelled
# hypothesis-generating because that is all a screen can be. The confirmatory
# question, which terms explain the data when they compete in one model, is
# ap_explains(), and the two are kept apart on purpose.

#' Screen every metadata variable against alpha and beta diversity
#'
#' Tests each usable metadata variable against every alpha diversity metric and
#' every beta diversity distance, ranks the results by effect size, and reports
#' how stable each rank is when the samples are subsampled.
#'
#' @section Ranking:
#' Rows are ordered by variance explained adjusted for degrees of freedom,
#' never by p. Raw R2 is biased upward: with no effect at all its expectation is
#' about `df / (n - 1)`, so a variable with many levels, or one recorded on few
#' samples, outranks a real two-group effect on the raw scale. Beta rows are
#' therefore ranked on adjusted PERMANOVA R2, `1 - (1 - R2)(n - 1)/(n - df - 1)`,
#' and numeric alpha rows on the same adjustment of Spearman rho squared.
#' Categorical alpha rows use epsilon squared (Kruskal-Wallis), which already
#' subtracts the null expectation of H. Raw values stay in `effect`, adjusted
#' ones in `effect_adj`.
#'
#' These are all proportions of variance, but not of the same variance: R2
#' partitions distance sums of squares and epsilon squared partitions rank
#' variance. Order within a family is exact; order across families is
#' approximate.
#'
#' @section Multiple testing:
#' One Benjamini-Hochberg correction is applied across every test the screen
#' ran, alpha and beta together. The number of tests is part of the result and
#' is printed with it.
#'
#' @section Stability:
#' `n_resample` subsamples of `fraction` of the samples are drawn without
#' replacement, and every effect size is recomputed on each. Stability is the
#' share of subsamples in which a row ranks in the top `top_k`. Subsampling is
#' used rather than the bootstrap because resampling with replacement puts
#' duplicate samples at zero distance from each other, which inflates
#' PERMANOVA R2. Subsamples overlap, so stability says how much a rank depends
#' on which samples were drawn. It does not say whether the effect is real.
#'
#' @section Variables that are not screened:
#' Constant variables, identifiers, variables with more than `max_levels`
#' levels, repeated-measures identifiers (which need a mixed model, not a
#' screen), and categorical variables whose smallest level has fewer than 3
#' samples. Each is listed in `skipped` with the reason.
#'
#' @param alpha Optional `ap_alpha` object from [ap_alpha()].
#' @param beta Optional `ap_beta` object from [ap_beta()]. At least one of
#'   `alpha` and `beta` is required.
#' @param variables Metadata variables to screen. `NULL` (default) screens every
#'   usable variable. Named variables still pass through the same guards.
#' @param permutations Permutations for each PERMANOVA and dispersion test.
#'   Default `999`.
#' @param n_resample Subsamples for the stability estimate. Default `100`,
#'   minimum `10`.
#' @param fraction Share of samples in each subsample, strictly between 0 and
#'   1. Default `0.8`.
#' @param top_k A row is stable in a subsample when it ranks in the top `top_k`.
#'   Default `5`.
#' @param seed Random seed for the permutations and the subsamples, recorded
#'   with the result.
#' @param max_levels Passed to [ap_scan_metadata()]. Default `20`.
#'
#' @return An object of class `ap_screen`: a list with `results` (one row per
#'   test, ranked by effect), `by_variable`, `skipped`, `untestable`,
#'   `n_tests`, `expected_false_positives`, `resamples`, and the settings used.
#' @export
ap_screen <- function(alpha = NULL,
                      beta = NULL,
                      variables = NULL,
                      permutations = 999L,
                      n_resample = 100L,
                      fraction = 0.8,
                      top_k = 5L,
                      seed = 1L,
                      max_levels = 20L) {
  ap_assert(!is.null(alpha) || !is.null(beta),
            "Give at least one of `alpha` (from `ap_alpha()`) or `beta` (from `ap_beta()`).")
  if (!is.null(alpha)) {
    ap_assert(inherits(alpha, "ap_alpha"),
              "`alpha` must come from `ap_alpha()`, not {class(alpha)[1]}.")
  }
  if (!is.null(beta)) {
    ap_assert(inherits(beta, "ap_beta"),
              "`beta` must come from `ap_beta()`, not {class(beta)[1]}.")
  }
  ap_assert(is.numeric(fraction) && length(fraction) == 1L && fraction > 0 && fraction < 1,
            paste0("`fraction` must be strictly between 0 and 1, not {fraction}. At 1 every ",
                   "subsample is the full data and stability is 1 by construction."))
  ap_assert(n_resample >= 10L,
            paste0("`n_resample` must be at least 10, not {n_resample}. Stability is part of ",
                   "the result, and fewer subsamples cannot estimate it."))

  meta <- if (!is.null(alpha)) alpha$metadata else beta$metadata
  sel <- ap_screen_variables(meta, variables, max_levels)
  ap_assert(nrow(sel$keep) > 0L,
            "No variable is testable. Skipped: {paste(sel$skipped$variable, collapse = ', ')}.")

  # Categorical variables become factors once, here, so the alpha tests, the
  # PERMANOVA and the subsamples all see the same coding. A numeric 0/1 column
  # is a grouping, not a gradient.
  coded <- lapply(seq_len(nrow(sel$keep)), function(i) {
    v <- meta[[sel$keep$variable[i]]]
    if (sel$keep$type[i] == "categorical") factor(as.character(v)) else as.numeric(v)
  })
  names(coded) <- sel$keep$variable

  rows <- list()
  specs <- list()
  untestable <- list()

  if (!is.null(alpha)) {
    for (m in unique(alpha$values$metric)) {
      vv <- alpha$values[alpha$values$metric == m, ]
      y <- stats::setNames(vv$value, vv$sample_id)
      for (i in seq_len(nrow(sel$keep))) {
        v <- sel$keep$variable[i]
        type <- sel$keep$type[i]
        g <- coded[[v]][match(names(y), rownames(meta))]
        e <- ap_screen_alpha_effect(y, g, type)
        if (is.na(e$p)) {
          untestable[[length(untestable) + 1L]] <- data.frame(
            family = "alpha", metric = m, variable = v,
            reason = "too few non-missing values or no variation", stringsAsFactors = FALSE)
          next
        }
        rows[[length(rows) + 1L]] <- data.frame(
          family = "alpha", metric = m, variable = v, type = type, test = e$test,
          n = e$n, df = e$df, statistic = e$statistic, effect_name = e$effect_name,
          effect = e$effect, effect_adj = e$effect_adj, p = e$p,
          dispersion_p = NA_real_, verdict = NA_character_,
          stringsAsFactors = FALSE)
        specs[[length(specs) + 1L]] <- list(family = "alpha", metric = m, ids = names(y),
                                            y = unname(y), g = g, type = type)
      }
    }
  }

  mats <- NULL
  if (!is.null(beta)) {
    mats <- lapply(beta$distances, as.matrix)
    beta_local <- beta
    bmeta <- beta$metadata
    for (v in names(coded)) {
      beta_local$metadata[[v]] <- coded[[v]][match(rownames(bmeta), rownames(meta))]
    }
    for (i in seq_len(nrow(sel$keep))) {
      v <- sel$keep$variable[i]
      type <- sel$keep$type[i]
      pn <- tryCatch(
        suppressMessages(ap_permanova(beta_local, v, permutations = permutations,
                                      by = "margin", seed = seed)),
        error = function(e) e)
      if (inherits(pn, "error")) {
        untestable[[length(untestable) + 1L]] <- data.frame(
          family = "beta", metric = paste(beta$metrics, collapse = ", "), variable = v,
          reason = conditionMessage(pn), stringsAsFactors = FALSE)
        next
      }
      for (m in beta$metrics) {
        r <- pn$results[pn$results$metric == m, ]
        d <- pn$dispersion[pn$dispersion$metric == m, ]
        it <- pn$interpretation[pn$interpretation$metric == m, ]
        ids <- rownames(mats[[m]])
        g <- beta_local$metadata[[v]][match(ids, rownames(bmeta))]

        # The subsamples use a direct trace formula instead of adonis2. If the
        # two ever disagree on the full data, every stability number is wrong
        # while still looking plausible, so the screen stops here.
        r2_fast <- ap_screen_r2(g, m = mats[[m]])
        ap_assert(
          isTRUE(abs(r2_fast - r$R2) < 1e-8),
          paste0("Internal R2 ({signif(r2_fast, 8)}) disagrees with adonis2 ",
                 "({signif(r$R2, 8)}) for `{v}` on {m}. Stability would be computed ",
                 "on the wrong quantity, so the screen stops.")
        )

        # The subsamples derive df from the grouping; the full-data df comes from
        # adonis2. They must agree or the adjusted scale differs between the two.
        ap_assert(isTRUE(r$df == ap_screen_df(g)),
                  "adonis2 reports {r$df} df for `{v}` on {m} but the grouping implies {ap_screen_df(g)}.")

        rows[[length(rows) + 1L]] <- data.frame(
          family = "beta", metric = m, variable = v, type = type,
          test = "PERMANOVA", n = r$n, df = r$df, statistic = r$pseudo_F, effect_name = "R2",
          effect = r$R2, effect_adj = ap_adjust_r2(r$R2, r$n, r$df), p = r$p,
          dispersion_p = if (nrow(d) == 0L) NA_real_ else d$dispersion_p[1],
          verdict = if (nrow(it) == 0L) NA_character_ else it$verdict[1],
          stringsAsFactors = FALSE)
        specs[[length(specs) + 1L]] <- list(family = "beta", metric = m, ids = ids,
                                            g = g, type = type)
      }
    }
  }

  ap_assert(length(rows) > 0L, "No test could be run on these variables.")
  results <- do.call(rbind, rows)
  results$q <- stats::p.adjust(results$p, method = "BH")
  n_tests <- nrow(results)
  ap_assert(top_k < n_tests,
            paste0("`top_k` ({top_k}) must be smaller than the number of tests ({n_tests}). ",
                   "Otherwise every test is in the top k and stability is 1 by construction."))

  ord <- order(results$effect_adj, decreasing = TRUE)
  results <- results[ord, , drop = FALSE]
  specs <- specs[ord]

  stab <- ap_screen_stability(specs, mats, n_resample, fraction, top_k, seed)
  results$stability <- stab$stability
  results$median_rank <- stab$median_rank
  results$rank <- seq_len(n_tests)
  rownames(results) <- NULL
  results <- results[, c("rank", "family", "metric", "variable", "type", "test", "n", "df",
                         "statistic", "effect_name", "effect", "effect_adj", "p", "q", "stability",
                         "median_rank", "dispersion_p", "verdict")]

  by_variable <- do.call(rbind, lapply(unique(results$variable), function(v) {
    j <- which(results$variable == v)
    in_top <- !is.na(stab$ranks[j, , drop = FALSE]) & stab$ranks[j, , drop = FALSE] <= top_k
    data.frame(
      variable = v, type = results$type[j[1]],
      best_family = results$family[j[1]], best_metric = results$metric[j[1]],
      best_effect_name = results$effect_name[j[1]], best_effect = results$effect_adj[j[1]],
      best_effect_raw = results$effect[j[1]],
      n_tests = length(j), n_q05 = sum(results$q[j] < 0.05),
      stability = mean(colSums(in_top) > 0),
      stringsAsFactors = FALSE)
  }))

  structure(
    list(results = results,
         by_variable = by_variable,
         skipped = sel$skipped,
         untestable = if (length(untestable)) do.call(rbind, untestable) else NULL,
         n_tests = n_tests,
         n_alpha = sum(results$family == "alpha"),
         n_beta = sum(results$family == "beta"),
         expected_false_positives = 0.05 * n_tests,
         resamples = stab$resamples,
         permutations = permutations, n_resample = n_resample, fraction = fraction,
         top_k = top_k, seed = seed, max_levels = max_levels),
    class = "ap_screen"
  )
}

#' @keywords internal
ap_screen_variables <- function(meta, variables, max_levels) {
  if (!is.null(variables)) {
    bad <- setdiff(variables, names(meta))
    ap_assert(length(bad) == 0L,
              "{cli::qty(length(bad))}Variable{?s} not in the sample metadata: {paste(bad, collapse = ', ')}.")
  }
  scan <- ap_scan_metadata(meta, max_levels = max_levels)
  cand <- variables %||% scan$variables$variable

  keep <- list()
  skipped <- list()
  for (v in cand) {
    r <- scan$variables[scan$variables$variable == v, ]
    reason <- if (v %in% scan$repeated_measures) {
      "values repeat across samples like a subject ID; samples are not independent, use a mixed model"
    } else if (r$role == "constant") {
      "constant, nothing to test"
    } else if (r$role == "identifier") {
      "one value per sample, an identifier"
    } else if (r$role == "high-cardinality") {
      sprintf("more than %d levels", max_levels)
    } else if (r$role %in% c("binary", "categorical") && !is.na(r$smallest_level_n) &&
               r$smallest_level_n < 3L) {
      sprintf("smallest level has %d sample(s); fewer than 3 cannot carry a test",
              r$smallest_level_n)
    } else {
      NA_character_
    }

    if (is.na(reason)) {
      type <- if (r$role == "numeric") "numeric" else "categorical"
      keep[[length(keep) + 1L]] <- data.frame(variable = v, type = type,
                                              stringsAsFactors = FALSE)
    } else {
      skipped[[length(skipped) + 1L]] <- data.frame(variable = v, reason = reason,
                                                    stringsAsFactors = FALSE)
    }
  }

  list(
    keep = if (length(keep)) do.call(rbind, keep) else
      data.frame(variable = character(0), type = character(0), stringsAsFactors = FALSE),
    skipped = if (length(skipped)) do.call(rbind, skipped) else NULL
  )
}

#' @keywords internal
ap_screen_alpha_effect <- function(y, g, type) {
  keep <- !is.na(y) & !is.na(g)
  y <- y[keep]
  g <- g[keep]
  n <- length(y)
  categorical <- type == "categorical"
  out <- list(test = if (categorical) "Kruskal-Wallis" else "Spearman",
              effect_name = if (categorical) "epsilon squared" else "rho squared",
              n = n, df = NA_real_, statistic = NA_real_, effect = NA_real_,
              effect_adj = NA_real_, p = NA_real_)

  if (categorical) {
    g <- droplevels(factor(g))
    k <- nlevels(g)
    if (k < 2L || n <= k || length(unique(y)) < 2L) return(out)
    ht <- tryCatch(stats::kruskal.test(y, g), error = function(e) NULL)
    if (is.null(ht)) return(out)
    out$statistic <- unname(ht$statistic)
    out$effect <- ap_epsilon_squared_point(out$statistic, n, k)
    out$df <- k - 1
    # Epsilon squared already subtracts the null expectation of H, which is
    # k - 1, so it needs no further adjustment.
    out$effect_adj <- out$effect
    out$p <- ht$p.value
  } else {
    if (n < 4L || length(unique(y)) < 2L || length(unique(g)) < 2L) return(out)
    ct <- tryCatch(suppressWarnings(stats::cor.test(y, g, method = "spearman", exact = FALSE)),
                   error = function(e) NULL)
    if (is.null(ct)) return(out)
    out$statistic <- unname(ct$estimate)
    out$effect <- out$statistic^2
    out$df <- 1
    out$effect_adj <- ap_adjust_r2(out$effect, n, 1)
    out$p <- ct$p.value
  }
  out
}

# Gower's double-centred matrix, G = -1/2 J D^2 J. PERMANOVA is built on it:
# total sum of squares is tr(G), and the sum of squares a model explains is
# tr(HG), with H the hat matrix of the model. This is the McArdle and Anderson
# formulation adonis2 uses.
#' @keywords internal
ap_gower <- function(m) {
  a <- -0.5 * m^2
  rm <- rowMeans(a)
  a - outer(rm, rm, "+") + mean(a)
}

#' @keywords internal
ap_model_r2 <- function(G, X) {
  qx <- qr(X)
  q <- qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]
  sum((G %*% q) * q) / sum(diag(G))
}

# R2 of one term on a distance matrix. `G` may be passed in precomputed for
# speed, but it is only valid for the full sample set it was built on: when `g`
# has missing values the centring changes, so G is rebuilt on the complete
# samples from `m`. Reusing it there would give a wrong number that looks right.
#' @keywords internal
ap_screen_r2 <- function(g, G = NULL, m = NULL) {
  keep <- !is.na(g)
  if (!all(keep) || is.null(G)) {
    ap_assert(!is.null(m), "A distance matrix is needed to rebuild the Gower matrix.")
    m <- m[keep, keep, drop = FALSE]
    g <- g[keep]
    G <- ap_gower(m)
  }
  if (is.factor(g)) {
    g <- droplevels(g)
    if (nlevels(g) < 2L) return(NA_real_)
  } else if (length(unique(g)) < 2L) {
    return(NA_real_)
  }
  ap_model_r2(G, stats::model.matrix(~ g))
}

# Adjusted R2. Raw R2 has a null expectation of about df / (n - 1), which is
# what lets a many-level or small-n variable outrank a real effect.
#' @keywords internal
ap_adjust_r2 <- function(r2, n, df) {
  denom <- n - df - 1
  ifelse(denom > 0, 1 - (1 - r2) * (n - 1) / denom, NA_real_)
}

# Degrees of freedom a single term uses, from its non-missing values.
#' @keywords internal
ap_screen_df <- function(g) {
  g <- g[!is.na(g)]
  if (is.factor(g)) nlevels(droplevels(g)) - 1L else 1L
}

#' @keywords internal
ap_screen_stability <- function(specs, mats, n_resample, fraction, top_k, seed) {
  universe <- sort(unique(unlist(lapply(specs, `[[`, "ids"))))
  size <- floor(fraction * length(universe))
  set.seed(seed)
  resamples <- lapply(seq_len(n_resample), function(b) {
    sample(universe, size, replace = FALSE)
  })

  n_tests <- length(specs)
  eff <- matrix(NA_real_, nrow = n_tests, ncol = n_resample)
  for (b in seq_len(n_resample)) {
    sub <- resamples[[b]]
    gower_cache <- list()
    for (j in seq_len(n_tests)) {
      s <- specs[[j]]
      idx <- which(s$ids %in% sub)
      if (s$family == "alpha") {
        eff[j, b] <- ap_screen_alpha_effect(s$y[idx], s$g[idx], s$type)$effect_adj
      } else {
        msub <- mats[[s$metric]][idx, idx, drop = FALSE]
        if (is.null(gower_cache[[s$metric]])) gower_cache[[s$metric]] <- ap_gower(msub)
        gs <- s$g[idx]
        r2 <- ap_screen_r2(gs, G = gower_cache[[s$metric]], m = msub)
        eff[j, b] <- ap_adjust_r2(r2, sum(!is.na(gs)), ap_screen_df(gs))
      }
    }
  }

  ranks <- apply(eff, 2, function(e) rank(-e, ties.method = "min", na.last = "keep"))
  ranks <- matrix(ranks, nrow = n_tests)
  list(
    stability = rowMeans(!is.na(ranks) & ranks <= top_k),
    median_rank = apply(ranks, 1, stats::median, na.rm = TRUE),
    ranks = ranks,
    resamples = resamples
  )
}

#' @keywords internal
ap_effect_abbrev <- function(name) {
  unname(c("R2" = "R2", "epsilon squared" = "eps2", "rho squared" = "rho2")[name])
}

#' @export
print.ap_screen <- function(x, n = 20L, ...) {
  cli::cli_h1("Exploratory screen (hypothesis-generating)")
  msg <- sprintf("%d tests (%d alpha, %d beta) across %d variables. One BH correction across all %d.",
                 x$n_tests, x$n_alpha, x$n_beta, nrow(x$by_variable), x$n_tests)
  cli::cli_text("{msg}")
  msg <- sprintf(paste0("If nothing were real, about %.0f of these would still reach p < 0.05. ",
                        "Read q, and read the effect size before either."),
                 x$expected_false_positives)
  cli::cli_alert_warning("{msg}")
  msg <- sprintf(paste0("Ranked by variance explained adjusted for degrees of freedom, not by p. ",
                        "Stability is the share of %d ",
                        "subsamples (%.0f%% of samples, without replacement, seed %s) in which ",
                        "the row ranks in the top %d."),
                 x$n_resample, 100 * x$fraction, x$seed, x$top_k)
  cli::cli_text("{msg}")

  r <- utils::head(x$results, n)
  out <- data.frame(
    rank = r$rank,
    variable = r$variable,
    family = r$family,
    metric = r$metric,
    adjusted = sprintf("%s %.4f", ap_effect_abbrev(r$effect_name), r$effect_adj),
    raw = sprintf("%.4f", r$effect),
    p = format.pval(r$p, digits = 2),
    q = format.pval(r$q, digits = 2),
    stable = sprintf("%.0f%%", 100 * r$stability),
    verdict = ifelse(is.na(r$verdict), "-", r$verdict),
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE)
  if (nrow(x$results) > n) {
    more <- nrow(x$results) - n
    cli::cli_text("{more} more row{?s} in {.code $results}.")
  }

  cli::cli_h2("By variable")
  bv <- x$by_variable
  print(data.frame(
    variable = bv$variable,
    best = sprintf("%s %.4f (%s, %s)", ap_effect_abbrev(bv$best_effect_name), bv$best_effect,
                   bv$best_family, bv$best_metric),
    `q<0.05` = sprintf("%d/%d", bv$n_q05, bv$n_tests),
    stable = sprintf("%.0f%%", 100 * bv$stability),
    check.names = FALSE, stringsAsFactors = FALSE
  ), row.names = FALSE)

  if (!is.null(x$skipped)) {
    cli::cli_h2("Not screened")
    for (i in seq_len(nrow(x$skipped))) {
      cli::cli_li("{.field {x$skipped$variable[i]}}: {x$skipped$reason[i]}")
    }
  }
  if (!is.null(x$untestable)) {
    cli::cli_h2("Could not be tested")
    for (i in seq_len(nrow(x$untestable))) {
      u <- x$untestable[i, ]
      cli::cli_li("{.field {u$variable}} ({u$family}, {u$metric}): {u$reason}")
    }
  }

  cli::cli_text("")
  cli::cli_alert_info(paste0(
    "Raw R2 and rho2 sit near df/(n-1) even with no effect, which favours many-level and ",
    "small-n variables. The ranking uses the adjusted values; raw ones are shown for reference."
  ))
  cli::cli_alert_info(paste0(
    "R2 partitions distance sums of squares and eps2 or rho2 partition rank variance. ",
    "Both are variance explained, not the same variance, so order across families is approximate."
  ))
  cli::cli_alert_info(paste0(
    "Subsamples overlap by design. Stability shows how much a rank depends on which samples ",
    "were drawn; it does not show that an effect is real."
  ))
  cli::cli_text(
    "To ask which variables explain the data when they compete in one model, use {.fn ap_explains}."
  )
  invisible(x)
}
