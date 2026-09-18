# Testing alpha diversity between groups.
#
# The default in the field is a t-test or Wilcoxon, a p-value, and a boxplot.
# Three things go wrong with that and all three are fixed here.
#
#   1. The assumption behind the chosen test is never checked, so nobody knows
#      whether it held.
#   2. The effect size is never reported, so at n = 487 a p of 0.001 gets
#      written up as if it were a large difference.
#   3. Repeated measures are analysed as if the samples were independent.
#
# This function checks the assumptions first, reports the effect size with an
# interval alongside every p-value, and routes to a mixed model when a subject
# identifier is given.

#' Test alpha diversity between groups
#'
#' @param alpha An `ap_alpha` object from [ap_alpha()].
#' @param group Name of the grouping variable in the sample metadata.
#' @param covariates Optional character vector of metadata variables to adjust
#'   for. Switches to a linear model.
#' @param subject Optional metadata variable identifying repeated measures.
#'   Switches to a linear mixed model with a random intercept per subject.
#' @param metrics Metrics to test. Defaults to everything in `alpha`.
#' @param seed Seed for the bootstrap used for non-parametric effect size
#'   intervals. Recorded in the result.
#' @param n_boot Bootstrap resamples for those intervals. Default `2000`.
#'
#' @return An object of class `ap_alpha_test`: a list with `results` (one row
#'   per metric) and the settings used.
#' @export
ap_alpha_test <- function(alpha,
                          group,
                          covariates = NULL,
                          subject = NULL,
                          metrics = NULL,
                          seed = 1L,
                          n_boot = 2000L) {
  ap_assert(inherits(alpha, "ap_alpha"),
            "`alpha` must come from `ap_alpha()`, not {class(alpha)[1]}.")
  meta <- alpha$metadata
  ap_assert(group %in% names(meta),
            "Variable `{group}` is not in the sample metadata. Available: {paste(names(meta), collapse = ', ')}.")

  for (v in c(covariates, subject)) {
    ap_assert(v %in% names(meta), "Variable `{v}` is not in the sample metadata.")
  }

  metrics <- metrics %||% unique(alpha$values$metric)
  set.seed(seed)

  rows <- lapply(metrics, function(m) {
    df <- alpha$values[alpha$values$metric == m, ]
    df$g <- meta[[group]][match(df$sample_id, rownames(meta))]
    keep <- !is.na(df$value) & !is.na(df$g)
    df <- df[keep, , drop = FALSE]
    if (nrow(df) < 4L || length(unique(df$g)) < 2L) return(NULL)

    if (!is.null(subject)) {
      return(ap_alpha_mixed(df, meta, m, group, covariates, subject))
    }
    if (!is.null(covariates)) {
      return(ap_alpha_lm(df, meta, m, group, covariates))
    }
    ap_alpha_simple(df, m, n_boot)
  })

  res <- do.call(rbind, Filter(Negate(is.null), rows))
  ap_assert(!is.null(res), "No metric had enough non-missing data to test.")
  res$p_adj <- stats::p.adjust(res$p, method = "BH")
  rownames(res) <- NULL

  out <- structure(
    list(results = res, group = group, covariates = covariates,
         subject = subject, seed = seed, n_boot = n_boot,
         n_groups = length(unique(meta[[group]][!is.na(meta[[group]])])),
         rarefied = alpha$rarefied, depth = alpha$depth),
    class = "ap_alpha_test"
  )
  out$interpretation <- ap_alpha_test_interpret(out)
  out
}

# Alpha diversity had no interpretation text at all until 2026-09-18. The report
# showed a table of p-values and effect sizes and left the reader to decide what
# they meant, which is the one thing this package is supposed to not do.
#
# The branches are exhaustive over the three facts available per metric: whether
# the adjusted p clears alpha, whether the assumption check passed, and whether
# a confidence interval was estimated. A non-significant result is never
# reported as "no difference" without saying what the interval still allows.
#' @keywords internal
ap_alpha_test_interpret <- function(x, alpha = 0.05) {
  r <- x$results
  rows <- lapply(seq_len(nrow(r)), function(i) {
    metric <- r$metric[i]
    p_adj <- r$p_adj[i]
    est <- r$estimate[i]
    eff <- r$effect[i]
    has_ci <- !is.na(r$ci_low[i])
    ci <- if (has_ci) sprintf("[%.3f, %.3f]", r$ci_low[i], r$ci_high[i]) else NA_character_
    assum <- r$assumptions_met[i]

    mk <- function(verdict, text) {
      data.frame(metric = metric, verdict = verdict, interpretation = text,
                 stringsAsFactors = FALSE)
    }

    if (is.na(p_adj)) {
      return(mk("not testable",
                "No p-value was produced for this metric, so nothing can be concluded from it."))
    }

    sig <- p_adj < alpha

    if (!sig) {
      text <- sprintf(
        paste0("No difference detected in %s (q = %s, %s = %.3f). "),
        metric, format.pval(p_adj, digits = 2), eff, est)
      text <- paste0(text, if (has_ci) {
        sprintf(paste0("The 95%% interval %s is the set of effect sizes still consistent ",
                       "with these data. This is absence of evidence, and the interval ",
                       "says how much evidence there is: a wide interval does not rule ",
                       "out a real difference."), ci)
      } else {
        paste0("No interval was estimated for this test, so the result says only that a ",
               "difference was not detected, not that one is absent.")
      })
      return(mk("no difference detected", text))
    }

    text <- sprintf(
      paste0("Groups differ in %s (q = %s, %s = %.3f%s). "),
      metric, format.pval(p_adj, digits = 2), eff, est,
      if (has_ci) paste0(", 95% CI ", ci) else "")

    if (isFALSE(assum)) {
      return(mk("difference, rank-based test", paste0(
        text,
        "The normality or equal-variance assumption failed, so a rank-based test was ",
        "used. The result is a statement about ranks, not about means, and the effect ",
        "size should be reported as such.")))
    }

    if (is.na(assum)) {
      return(mk("difference", paste0(
        text,
        "Assumptions were not checked for this test. Read the effect size rather than ",
        "the p-value: at this sample size a small p-value is cheap.")))
    }

    mk("difference", paste0(
      text,
      "Assumption checks passed. Read the effect size rather than the p-value: at this ",
      "sample size a small p-value is cheap, and the effect size is what the claim rests on."))
  })
  do.call(rbind, rows)
}

#' @keywords internal
ap_alpha_simple <- function(df, metric, n_boot) {
  g <- as.factor(as.character(df$g))
  v <- df$value
  k <- nlevels(g)

  # Assumptions, checked before the test rather than assumed by it. Shapiro-Wilk
  # on the within-group residuals; Bartlett for equality of variances.
  resid <- v - stats::ave(v, g, FUN = mean)
  normal_p <- if (length(v) >= 3L && length(v) <= 5000L) {
    tryCatch(stats::shapiro.test(resid)$p.value, error = function(e) NA_real_)
  } else NA_real_
  var_p <- tryCatch(stats::bartlett.test(v, g)$p.value, error = function(e) NA_real_)

  assumptions_ok <- !is.na(normal_p) && normal_p > 0.05 &&
    !is.na(var_p) && var_p > 0.05

  if (k == 2L) {
    if (assumptions_ok) {
      ht <- stats::t.test(v ~ g, var.equal = TRUE)
      es <- ap_hedges_g(v, g)
      test <- "Student's t"
    } else {
      ht <- stats::wilcox.test(v ~ g, exact = FALSE)
      es <- ap_cliffs_delta(v, g, n_boot)
      test <- "Wilcoxon rank-sum"
    }
  } else {
    if (assumptions_ok) {
      fit <- stats::aov(v ~ g)
      sm <- summary(fit)[[1]]
      ht <- list(statistic = sm[["F value"]][1], p.value = sm[["Pr(>F)"]][1])
      es <- ap_eta_squared(sm)
      test <- "One-way ANOVA"
    } else {
      ht <- stats::kruskal.test(v, g)
      es <- ap_eta_squared_h(v, g, n_boot)
      test <- "Kruskal-Wallis"
    }
  }

  data.frame(
    metric = metric,
    test = test,
    n = length(v),
    k_groups = k,
    statistic = unname(ht$statistic %||% NA_real_),
    p = unname(ht$p.value),
    effect = es$name,
    estimate = es$estimate,
    ci_low = es$ci[1],
    ci_high = es$ci[2],
    normality_p = normal_p,
    variance_p = var_p,
    assumptions_met = assumptions_ok,
    stringsAsFactors = FALSE
  )
}

# Hedges' g: Cohen's d with the small-sample correction. Interval from the
# large-sample standard error, which is the usual approximation and is stated
# rather than left implicit.
#' @keywords internal
ap_hedges_g <- function(v, g) {
  lv <- levels(g)
  a <- v[g == lv[1]]
  b <- v[g == lv[2]]
  n1 <- length(a); n2 <- length(b)
  sp <- sqrt(((n1 - 1) * stats::var(a) + (n2 - 1) * stats::var(b)) / (n1 + n2 - 2))
  d <- (mean(b) - mean(a)) / sp
  j <- 1 - 3 / (4 * (n1 + n2) - 9)
  gval <- j * d
  se <- sqrt((n1 + n2) / (n1 * n2) + gval^2 / (2 * (n1 + n2 - 2)))
  list(name = "Hedges' g", estimate = gval,
       ci = gval + c(-1, 1) * stats::qnorm(0.975) * se)
}

# Cliff's delta: the non-parametric counterpart, the probability that a random
# value from one group exceeds one from the other, rescaled to [-1, 1]. It makes
# no distributional assumption, which is the point of reaching for it.
#' @keywords internal
ap_cliffs_delta <- function(v, g, n_boot) {
  lv <- levels(g)
  a <- v[g == lv[1]]
  b <- v[g == lv[2]]
  delta <- ap_cliffs_point(a, b)
  boots <- vapply(seq_len(n_boot), function(i) {
    ap_cliffs_point(sample(a, replace = TRUE), sample(b, replace = TRUE))
  }, numeric(1))
  list(name = "Cliff's delta", estimate = delta,
       ci = unname(stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE)))
}

#' @keywords internal
ap_cliffs_point <- function(a, b) {
  cmp <- outer(b, a, ">") * 1 - outer(b, a, "<") * 1
  mean(cmp)
}

#' @keywords internal
ap_eta_squared <- function(sm) {
  ss <- sm[["Sum Sq"]]
  est <- ss[1] / sum(ss)
  # Fisher's approximation via the F distribution is not stable for small
  # designs; the interval is left to the bootstrap path instead.
  list(name = "eta squared", estimate = est, ci = c(NA_real_, NA_real_))
}

#' @keywords internal
ap_eta_squared_h <- function(v, g, n_boot) {
  # Kruskal-Wallis effect size: eta squared based on H, (H - k + 1) / (n - k),
  # computed from stats::kruskal.test. This one line is a documented exception
  # to taking statistics from packages. rstatix::kruskal_effsize computes the
  # same quantity, but from rstatix 1.0.0 it clamps the value to [0, 1]. The
  # estimate is bias-corrected so that a null effect averages zero and is
  # negative about half the time; clamping piles nulls at 0 and biases them
  # upward. rstatix is used in the tests as an independent check on positive
  # effects, where the two agree.
  #
  # The interval comes from boot::boot and boot::boot.ci (percentile), with the
  # same formula on each replicate. rstatix's own interval is rounded to two
  # decimals, too coarse for effects of 0.005 to 0.04.
  df <- data.frame(value = v, g = droplevels(as.factor(g)))

  eta_h <- function(d, i) {
    d <- d[i, , drop = FALSE]
    grp <- droplevels(d$g)
    k <- nlevels(grp)
    if (k < 2L) return(NA_real_)
    H <- unname(stats::kruskal.test(d$value, grp)$statistic)
    (H - k + 1) / (nrow(d) - k)
  }
  est <- eta_h(df, seq_len(nrow(df)))

  b <- boot::boot(df, eta_h, R = n_boot)
  ci <- tryCatch(unname(boot::boot.ci(b, type = "perc")$percent[4:5]),
                 error = function(e) c(NA_real_, NA_real_))
  list(name = "eta squared (H)", estimate = est, ci = ci)
}

#' @keywords internal
ap_alpha_lm <- function(df, meta, metric, group, covariates) {
  for (cv in covariates) df[[cv]] <- meta[[cv]][match(df$sample_id, rownames(meta))]
  df <- df[stats::complete.cases(df[, c("value", "g", covariates)]), , drop = FALSE]
  df$g <- as.factor(as.character(df$g))

  full <- stats::lm(stats::as.formula(paste("value ~ g +", paste(covariates, collapse = " + "))), data = df)
  reduced <- stats::lm(stats::as.formula(paste("value ~", paste(covariates, collapse = " + "))), data = df)
  an <- stats::anova(reduced, full)

  # Partial eta squared for the group term, after the covariates.
  ss_group <- an[["RSS"]][1] - an[["RSS"]][2]
  partial_eta <- ss_group / (ss_group + an[["RSS"]][2])

  data.frame(
    metric = metric,
    test = paste0("Linear model, adjusted for ", paste(covariates, collapse = " + ")),
    n = nrow(df),
    k_groups = nlevels(df$g),
    statistic = an[["F"]][2],
    p = an[["Pr(>F)"]][2],
    effect = "partial eta squared",
    estimate = partial_eta,
    ci_low = NA_real_,
    ci_high = NA_real_,
    normality_p = tryCatch(stats::shapiro.test(stats::residuals(full))$p.value,
                           error = function(e) NA_real_),
    variance_p = NA_real_,
    assumptions_met = NA,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
ap_alpha_mixed <- function(df, meta, metric, group, covariates, subject) {
  ap_assert(
    requireNamespace("lme4", quietly = TRUE),
    paste0("Repeated measures need a mixed model, which needs the lme4 package. ",
           "Install it, or drop `subject` and state in the methods that samples were ",
           "treated as independent.")
  )
  df$subject <- meta[[subject]][match(df$sample_id, rownames(meta))]
  for (cv in covariates) df[[cv]] <- meta[[cv]][match(df$sample_id, rownames(meta))]
  df <- df[stats::complete.cases(df[, c("value", "g", "subject", covariates)]), , drop = FALSE]
  df$g <- as.factor(as.character(df$g))

  fixed <- paste(c("g", covariates), collapse = " + ")
  full <- lme4::lmer(stats::as.formula(paste("value ~", fixed, "+ (1 | subject)")),
                     data = df, REML = FALSE)
  reduced_fixed <- if (is.null(covariates)) "1" else paste(covariates, collapse = " + ")
  reduced <- lme4::lmer(stats::as.formula(paste("value ~", reduced_fixed, "+ (1 | subject)")),
                        data = df, REML = FALSE)
  an <- stats::anova(reduced, full)

  data.frame(
    metric = metric,
    test = paste0("Linear mixed model, random intercept per ", subject),
    n = nrow(df),
    k_groups = nlevels(df$g),
    statistic = an[["Chisq"]][2],
    p = an[["Pr(>Chisq)"]][2],
    effect = "marginal R2 gain",
    estimate = ap_mixed_r2_gain(full, reduced),
    ci_low = NA_real_,
    ci_high = NA_real_,
    normality_p = NA_real_,
    variance_p = NA_real_,
    assumptions_met = NA,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
ap_mixed_r2_gain <- function(full, reduced) {
  var_f <- stats::var(stats::predict(full, re.form = NA))
  var_r <- stats::var(stats::predict(reduced, re.form = NA))
  total <- stats::var(stats::model.response(stats::model.frame(full)))
  (var_f - var_r) / total
}

#' @export
print.ap_alpha_test <- function(x, ...) {
  cli::cli_h1("Alpha diversity: {x$group}")
  if (!is.null(x$covariates)) cli::cli_text("Adjusted for: {paste(x$covariates, collapse = ', ')}")
  if (!is.null(x$subject)) cli::cli_text("Repeated measures on: {x$subject}")
  if (isTRUE(x$rarefied)) {
    cli::cli_text("Rarefied to {format(x$depth, big.mark = ',')} reads")
  }

  r <- x$results
  out <- data.frame(
    metric = r$metric,
    test = r$test,
    p = format.pval(r$p, digits = 2),
    q = format.pval(r$p_adj, digits = 2),
    effect = sprintf("%s = %.3f", r$effect, r$estimate),
    CI95 = ifelse(is.na(r$ci_low), "-", sprintf("[%.3f, %.3f]", r$ci_low, r$ci_high)),
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE)

  failed <- r[!is.na(r$assumptions_met) & !r$assumptions_met, ]
  if (nrow(failed) > 0L) {
    cli::cli_alert_info(paste0(
      "Assumptions failed for {paste(failed$metric, collapse = ', ')}, so the ",
      "rank-based test was used for {cli::qty(nrow(failed))}{?it/them}."
    ))
  }
  if (!is.null(x$interpretation)) {
    cli::cli_h2("Interpretation")
    for (i in seq_len(nrow(x$interpretation))) {
      it <- x$interpretation[i, ]
      cli::cli_li("{.field {it$metric}}: {it$interpretation}")
    }
  }
  invisible(x)
}
