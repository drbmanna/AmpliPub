# PERMANOVA, with the dispersion test bound in.
#
# PERMANOVA asks whether group centroids differ. It is sensitive to differences
# in within-group dispersion, so a significant result can mean the groups sit
# in different places, or that one group is simply more variable than the
# other, or both. `betadisper` distinguishes them, almost nobody runs it, and
# without it a significant PERMANOVA does not support the claim it is usually
# used to make.
#
# So the dispersion test is not an option here. It runs every time, and the
# result object cannot be printed without it. The interpretation is written out
# in words rather than left for the reader to reconstruct from two p-values.
#
# The other under-reported choice is `adonis2(by = )`. With `by = "terms"`,
# vegan attributes variance to terms sequentially, so reordering the formula
# changes every result but the last. AmpliPub defaults to `by = "margin"`,
# where each term is assessed after all the others and order does not matter,
# and records which was used.

#' PERMANOVA with a bound dispersion check
#'
#' Runs `vegan::adonis2` across every distance metric in an `ap_beta` object,
#' and `vegan::betadisper` with `vegan::permutest` for each categorical term,
#' then states what the combination means.
#'
#' @param beta An `ap_beta` object from [ap_beta()].
#' @param terms Character vector of metadata variables, or a one-sided formula
#'   such as `~ Site + dx`. Multiple terms are adjusted for one another.
#' @param metrics Metrics to test. Defaults to everything in `beta`.
#' @param permutations Number of permutations. Default `999`.
#' @param by `"margin"` (default) assesses each term after all others, so the
#'   formula order does not matter. `"terms"` is sequential and order-dependent.
#' @param strata Optional metadata variable to permute within, for a blocked
#'   design.
#' @param seed Random seed, recorded with the result.
#'
#' @return An object of class `ap_permanova`: a list with `results` (one row per
#'   metric per term), `dispersion` (one row per metric per categorical term),
#'   and the settings used.
#' @export
ap_permanova <- function(beta,
                         terms,
                         metrics = NULL,
                         permutations = 999L,
                         by = c("margin", "terms"),
                         strata = NULL,
                         seed = 1L) {
  ap_assert(inherits(beta, "ap_beta"),
            "`beta` must come from `ap_beta()`, not {class(beta)[1]}.")
  by <- match.arg(by)

  if (inherits(terms, "formula")) terms <- all.vars(terms)
  ap_assert(is.character(terms) && length(terms) > 0L,
            "`terms` must be variable names or a one-sided formula such as `~ Site + dx`.")

  meta <- beta$metadata
  for (v in c(terms, strata)) {
    ap_assert(v %in% names(meta),
              "Variable `{v}` is not in the sample metadata. Available: {paste(names(meta), collapse = ', ')}.")
  }

  metrics <- metrics %||% names(beta$distances)
  bad <- setdiff(metrics, names(beta$distances))
  ap_assert(length(bad) == 0L,
            "{cli::qty(length(bad))}Metric{?s} not in this ap_beta object: {paste(bad, collapse = ', ')}.")

  res_rows <- list()
  disp_rows <- list()

  for (m in metrics) {
    d <- beta$distances[[m]]
    ids <- attr(d, "Labels")
    md <- meta[ids, c(terms, strata), drop = FALSE]
    complete <- stats::complete.cases(md)

    if (sum(complete) < length(ids)) {
      cli::cli_inform(paste0(
        "{sum(!complete)} sample{?s} dropped from the {m} test for missing values in ",
        "{paste(c(terms, strata), collapse = ', ')}."
      ))
      d <- stats::as.dist(as.matrix(d)[complete, complete])
      md <- md[complete, , drop = FALSE]
    }
    # Characters become factors; numeric terms stay numeric and are fitted as
    # continuous predictors, which is what a numeric covariate should be.
    for (v in terms) if (!is.numeric(md[[v]])) md[[v]] <- factor(as.character(md[[v]]))

    set.seed(seed)
    fml <- stats::as.formula(paste("d ~", paste(terms, collapse = " + ")))
    control <- if (is.null(strata)) NULL else {
      permute::how(nperm = permutations, blocks = factor(as.character(md[[strata]])))
    }
    fit <- if (is.null(control)) {
      vegan::adonis2(fml, data = md, permutations = permutations, by = by)
    } else {
      vegan::adonis2(fml, data = md, permutations = control, by = by)
    }

    tbl <- as.data.frame(fit)
    keep <- rownames(tbl)[!rownames(tbl) %in% c("Residual", "Total")]
    for (tm in keep) {
      res_rows[[length(res_rows) + 1L]] <- data.frame(
        metric = m, term = tm, n = nrow(md),
        df = tbl[tm, "Df"],
        pseudo_F = tbl[tm, "F"],
        R2 = tbl[tm, "R2"],
        p = tbl[tm, "Pr(>F)"],
        stringsAsFactors = FALSE
      )
    }

    # Dispersion, per categorical term. betadisper needs a grouping factor, so a
    # continuous covariate has no dispersion test and is marked as such.
    for (v in terms) {
      if (is.numeric(md[[v]])) {
        disp_rows[[length(disp_rows) + 1L]] <- data.frame(
          metric = m, term = v, dispersion_F = NA_real_, dispersion_p = NA_real_,
          max_centroid_ratio = NA_real_,
          note = "continuous term; betadisper needs groups",
          stringsAsFactors = FALSE
        )
        next
      }
      g <- md[[v]]
      if (nlevels(droplevels(g)) < 2L) next
      set.seed(seed)
      bd <- vegan::betadisper(d, droplevels(g))
      pt <- vegan::permutest(bd, permutations = permutations)
      dist_means <- tapply(bd$distances, bd$group, mean)
      disp_rows[[length(disp_rows) + 1L]] <- data.frame(
        metric = m, term = v,
        dispersion_F = pt$tab[1, "F"],
        dispersion_p = pt$tab[1, "Pr(>F)"],
        max_centroid_ratio = max(dist_means) / min(dist_means),
        note = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }

  results <- do.call(rbind, res_rows)
  dispersion <- do.call(rbind, disp_rows)
  results$p_adj <- stats::p.adjust(results$p, method = "BH")
  rownames(results) <- NULL
  rownames(dispersion) <- NULL

  out <- structure(
    list(results = results, dispersion = dispersion, terms = terms,
         metrics = metrics, permutations = permutations, by = by,
         strata = strata, seed = seed, beta = beta),
    class = "ap_permanova"
  )
  out$interpretation <- ap_permanova_interpret(out)
  out
}

# The whole point of the function. Two p-values do not tell a reader what they
# are looking at unless someone says so.
#' @keywords internal
ap_permanova_interpret <- function(x, alpha = 0.05) {
  r <- x$results
  d <- x$dispersion
  rows <- lapply(seq_len(nrow(r)), function(i) {
    metric <- r$metric[i]
    term <- r$term[i]
    sig <- !is.na(r$p[i]) && r$p[i] < alpha
    hit <- if (is.null(d)) NULL else d[d$metric == metric & d$term == term, ]

    if (is.null(hit) || nrow(hit) == 0L || is.na(hit$dispersion_p[1])) {
      verdict <- if (sig) "location shift, dispersion not testable" else "no detectable difference"
      text <- if (sig) {
        "Significant, but this term is continuous, so betadisper cannot rule out a dispersion effect."
      } else {
        "No detectable difference in composition."
      }
      return(data.frame(metric = metric, term = term, verdict = verdict,
                        interpretation = text, stringsAsFactors = FALSE))
    }

    disp_sig <- hit$dispersion_p[1] < alpha
    ratio <- hit$max_centroid_ratio[1]

    if (sig && !disp_sig) {
      verdict <- "location shift"
      text <- sprintf(
        paste0("Groups differ in composition (pseudo-F = %.3f, R2 = %.4f, p = %s) and ",
               "dispersions do not differ (p = %s). This is a genuine shift in centroid, ",
               "which is what a PERMANOVA is normally taken to mean."),
        r$pseudo_F[i], r$R2[i], format.pval(r$p[i], digits = 2),
        format.pval(hit$dispersion_p[1], digits = 2))
    } else if (sig && disp_sig) {
      verdict <- "confounded by dispersion"
      text <- sprintf(
        paste0("PERMANOVA is significant (pseudo-F = %.3f, R2 = %.4f, p = %s) but so is ",
               "the dispersion test (p = %s, spread differs by %.2fx between groups). ",
               "PERMANOVA is sensitive to unequal dispersion, so this result cannot be ",
               "read as a centroid shift. It may be a shift, a spread difference, or both. ",
               "Report both numbers."),
        r$pseudo_F[i], r$R2[i], format.pval(r$p[i], digits = 2),
        format.pval(hit$dispersion_p[1], digits = 2), ratio)
    } else if (!sig && disp_sig) {
      verdict <- "dispersion only"
      text <- sprintf(
        paste0("Centroids do not differ (p = %s) but dispersions do (p = %s, %.2fx). ",
               "The groups sit in the same place and one is more variable. That is a ",
               "finding in its own right, not a null result."),
        format.pval(r$p[i], digits = 2), format.pval(hit$dispersion_p[1], digits = 2), ratio)
    } else {
      verdict <- "no difference"
      text <- sprintf(
        "Neither centroid (p = %s) nor dispersion (p = %s) differs.",
        format.pval(r$p[i], digits = 2), format.pval(hit$dispersion_p[1], digits = 2))
    }
    data.frame(metric = metric, term = term, verdict = verdict,
               interpretation = text, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' @export
print.ap_permanova <- function(x, ...) {
  cli::cli_h1("PERMANOVA with dispersion check")
  cli::cli_text("Terms: {paste(x$terms, collapse = ' + ')} | {x$permutations} permutations | ",
                "by = \"{x$by}\" | seed {x$seed}")
  if (!is.null(x$strata)) cli::cli_text("Permuted within: {x$strata}")
  if (isTRUE(x$beta$rarefied)) {
    cli::cli_text("Rarefied to {format(x$beta$depth, big.mark = ',')} reads")
  }

  tbl <- merge(x$results, x$dispersion[, c("metric", "term", "dispersion_F", "dispersion_p")],
               by = c("metric", "term"), all.x = TRUE, sort = FALSE)
  tbl <- merge(tbl, x$interpretation[, c("metric", "term", "verdict")],
               by = c("metric", "term"), all.x = TRUE, sort = FALSE)
  tbl <- tbl[order(match(tbl$metric, x$metrics), match(tbl$term, x$terms)), ]

  out <- data.frame(
    metric = tbl$metric,
    term = tbl$term,
    pseudo_F = round(tbl$pseudo_F, 3),
    R2 = round(tbl$R2, 4),
    p = format.pval(tbl$p, digits = 2),
    q = format.pval(tbl$p_adj, digits = 2),
    disp_p = ifelse(is.na(tbl$dispersion_p), "-", format.pval(tbl$dispersion_p, digits = 2)),
    verdict = tbl$verdict,
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE)

  cli::cli_h2("What this means")
  for (i in seq_len(nrow(x$interpretation))) {
    r <- x$interpretation[i, ]
    cli::cli_li("{.field {r$metric}} / {.field {r$term}}: {r$interpretation}")
  }

  cli::cli_text("")
  cli::cli_alert_info(paste0(
    "R2 is the number to report. At this sample size a small p-value costs almost ",
    "nothing; R2 says how much of the variation the term actually explains."
  ))
  invisible(x)
}
