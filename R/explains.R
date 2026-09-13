# The confirmatory question, kept apart from the screen.
#
# A screen tests each variable on its own, so two variables that carry the same
# signal both look strong. Put them in one model and each is assessed after the
# other: the variance only one of them can explain is its own, and the variance
# both could explain belongs to neither. On Baxter that is the Site and dx
# question. Reporting the shared part is the point; hiding it inside two
# marginal R2 values that do not add up to the model R2 is how an overlap gets
# written up as two findings.

#' Which terms explain the data when they compete in one model
#'
#' Fits every term together and reports, for each, the variance it explains
#' after all the others (marginal sums of squares), the variance explained by
#' the model as a whole, and the part of it no single term can claim.
#'
#' @section Shared variance:
#' `shared_R2 = model_R2 - sum(marginal R2)`. Positive values are variance the
#' terms explain jointly and cannot be assigned to any one of them, which is
#' what overlapping or confounded terms produce. Values near zero mean the terms
#' are close to independent. Negative values mean the terms suppress each other,
#' each explaining more once the other is in the model.
#'
#' @section Not a screen:
#' The terms here should be chosen before fitting. If they came out of
#' [ap_screen()], the methods section should say so.
#'
#' @param beta Optional `ap_beta` object from [ap_beta()].
#' @param alpha Optional `ap_alpha` object from [ap_alpha()]. At least one of
#'   `beta` and `alpha` is required.
#' @param terms Character vector of at least two metadata variables, or a
#'   one-sided formula such as `~ Site + dx`.
#' @param beta_metrics Distances to use. Defaults to everything in `beta`.
#' @param alpha_metrics Alpha metrics to use. Defaults to everything in `alpha`.
#' @param permutations Permutations for the PERMANOVA and dispersion tests.
#'   Default `999`.
#' @param seed Random seed, recorded with the result.
#'
#' @return An object of class `ap_explains`: a list with `beta` and `alpha`
#'   components, each holding `terms` (one row per metric per term) and `model`
#'   (one row per metric, with `model_R2`, `sum_marginal_R2` and `shared_R2`).
#' @export
ap_explains <- function(beta = NULL,
                        alpha = NULL,
                        terms,
                        beta_metrics = NULL,
                        alpha_metrics = NULL,
                        permutations = 999L,
                        seed = 1L) {
  ap_assert(!is.null(beta) || !is.null(alpha),
            "Give at least one of `beta` (from `ap_beta()`) or `alpha` (from `ap_alpha()`).")
  if (inherits(terms, "formula")) terms <- all.vars(terms)
  ap_assert(is.character(terms) && length(terms) >= 2L,
            paste0("`ap_explains()` weighs terms against each other, so it needs at least two. ",
                   "For a single variable use `ap_permanova()` or `ap_alpha_test()`."))

  out <- list(terms = terms, permutations = permutations, seed = seed)
  if (!is.null(beta)) {
    ap_assert(inherits(beta, "ap_beta"),
              "`beta` must come from `ap_beta()`, not {class(beta)[1]}.")
    out$beta <- ap_explains_beta(beta, terms, beta_metrics, permutations, seed)
  }
  if (!is.null(alpha)) {
    ap_assert(inherits(alpha, "ap_alpha"),
              "`alpha` must come from `ap_alpha()`, not {class(alpha)[1]}.")
    out$alpha <- ap_explains_alpha(alpha, terms, alpha_metrics)
  }
  structure(out, class = "ap_explains")
}

# Perfectly collinear terms have no marginal sum of squares to report: removing
# either leaves the fit unchanged, so both come back as zero and the model R2
# looks unexplained. That reads as a finding and is an artefact of the design.
#' @keywords internal
ap_check_collinear <- function(md, terms) {
  md <- md[stats::complete.cases(md[, terms, drop = FALSE]), terms, drop = FALSE]
  for (v in terms) if (!is.numeric(md[[v]])) md[[v]] <- factor(as.character(md[[v]]))
  X <- stats::model.matrix(~ ., data = md)
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    aliased <- colnames(X)[qx$pivot[(qx$rank + 1L):ncol(X)]]
    ap_abort(paste0(
      "These terms are perfectly collinear: {paste(aliased, collapse = ', ')} is fully ",
      "determined by the others. No term can be assessed after the others, so there is ",
      "nothing to partition. Drop one, and report that they cannot be separated in this design."
    ))
  }
  md
}

#' @keywords internal
ap_explains_beta <- function(beta, terms, metrics, permutations, seed) {
  meta <- beta$metadata
  for (v in terms) {
    ap_assert(v %in% names(meta),
              "Variable `{v}` is not in the sample metadata. Available: {paste(names(meta), collapse = ', ')}.")
  }
  ap_check_collinear(meta, terms)

  pn <- ap_permanova(beta, terms, metrics = metrics, permutations = permutations,
                     by = "margin", seed = seed)

  model <- do.call(rbind, lapply(pn$metrics, function(m) {
    d <- beta$distances[[m]]
    ids <- attr(d, "Labels")
    md <- meta[ids, terms, drop = FALSE]
    cc <- stats::complete.cases(md)
    md <- md[cc, , drop = FALSE]
    for (v in terms) if (!is.numeric(md[[v]])) md[[v]] <- factor(as.character(md[[v]]))
    mm <- as.matrix(d)[cc, cc, drop = FALSE]

    marg <- pn$results[pn$results$metric == m, ]
    ap_assert(all(marg$n == nrow(md)),
              "The model R2 and the marginal R2 for {m} were computed on different samples.")
    # R2 of all terms fitted together, from vegan::adonis2 with by = NULL.
    dd <- stats::as.dist(mm)
    model_r2 <- vegan::adonis2(dd ~ ., data = md, by = NULL, permutations = 0)[1, "R2"]
    data.frame(metric = m, n = nrow(md), model_R2 = model_r2,
               sum_marginal_R2 = sum(marg$R2), shared_R2 = model_r2 - sum(marg$R2),
               stringsAsFactors = FALSE)
  }))

  tbl <- merge(pn$results, pn$dispersion[, c("metric", "term", "dispersion_p")],
               by = c("metric", "term"), all.x = TRUE, sort = FALSE)
  tbl <- merge(tbl, pn$interpretation[, c("metric", "term", "verdict")],
               by = c("metric", "term"), all.x = TRUE, sort = FALSE)
  tbl <- tbl[order(match(tbl$metric, pn$metrics), match(tbl$term, terms)), ]
  rownames(tbl) <- NULL

  list(terms = tbl, model = model, permanova = pn)
}

#' @keywords internal
ap_explains_alpha <- function(alpha, terms, metrics) {
  meta <- alpha$metadata
  for (v in terms) {
    ap_assert(v %in% names(meta),
              "Variable `{v}` is not in the sample metadata. Available: {paste(names(meta), collapse = ', ')}.")
  }
  ap_check_collinear(meta, terms)
  metrics <- metrics %||% unique(alpha$values$metric)

  term_rows <- list()
  model_rows <- list()
  for (m in metrics) {
    vv <- alpha$values[alpha$values$metric == m, ]
    df <- data.frame(value = vv$value,
                     meta[match(vv$sample_id, rownames(meta)), terms, drop = FALSE],
                     check.names = FALSE, stringsAsFactors = FALSE)
    cc <- stats::complete.cases(df)
    if (any(!cc)) {
      n_drop <- sum(!cc)
      cli::cli_inform(paste0(
        "{n_drop} sample{?s} dropped from the {m} model for missing values in ",
        "{paste(terms, collapse = ', ')}."
      ))
    }
    df <- df[cc, , drop = FALSE]
    for (v in terms) if (!is.numeric(df[[v]])) df[[v]] <- factor(as.character(df[[v]]))
    ap_assert(nrow(df) > ncol(stats::model.matrix(~ ., data = df[, terms, drop = FALSE])) + 1L,
              "Too few complete samples ({nrow(df)}) to fit {length(terms)} terms for {m}.")

    fml <- stats::as.formula(paste("value ~", paste0("`", terms, "`", collapse = " + ")))
    fit <- stats::lm(fml, data = df)
    dr <- stats::drop1(fit, test = "F")
    sst <- sum((df$value - mean(df$value))^2)
    ss <- dr[["Sum of Sq"]][-1]
    model_r2 <- summary(fit)$r.squared

    term_rows[[length(term_rows) + 1L]] <- data.frame(
      metric = m, term = gsub("`", "", rownames(dr)[-1]), n = nrow(df),
      df = dr[["Df"]][-1], F = dr[["F value"]][-1], R2 = ss / sst, p = dr[["Pr(>F)"]][-1],
      stringsAsFactors = FALSE)
    res <- stats::residuals(fit)
    model_rows[[length(model_rows) + 1L]] <- data.frame(
      metric = m, n = nrow(df), model_R2 = model_r2, sum_marginal_R2 = sum(ss / sst),
      shared_R2 = model_r2 - sum(ss / sst),
      residual_normality_p = if (length(res) >= 3L && length(res) <= 5000L) {
        tryCatch(stats::shapiro.test(res)$p.value, error = function(e) NA_real_)
      } else NA_real_,
      stringsAsFactors = FALSE)
  }

  list(terms = do.call(rbind, term_rows), model = do.call(rbind, model_rows))
}

#' @export
print.ap_explains <- function(x, ...) {
  cli::cli_h1("Confirmatory model: {paste(x$terms, collapse = ' + ')}")
  cli::cli_text(paste0(
    "All terms in one model. Each term's R2 is what it explains after all the others ",
    "(marginal sums of squares), so term order does not matter."
  ))

  describe <- function(model) {
    for (i in seq_len(nrow(model))) {
      r <- model[i, ]
      msg <- sprintf("%s: model R2 %.4f; terms alone %.4f; shared %.4f",
                     r$metric, r$model_R2, r$sum_marginal_R2, r$shared_R2)
      cli::cli_li("{msg}")
    }
  }

  if (!is.null(x$beta)) {
    cli::cli_h2("Beta diversity (PERMANOVA, by = \"margin\", {x$permutations} permutations)")
    t <- x$beta$terms
    print(data.frame(
      metric = t$metric, term = t$term, R2 = round(t$R2, 4),
      pseudo_F = round(t$pseudo_F, 3), p = format.pval(t$p, digits = 2),
      disp_p = ifelse(is.na(t$dispersion_p), "-", format.pval(t$dispersion_p, digits = 2)),
      verdict = t$verdict, stringsAsFactors = FALSE
    ), row.names = FALSE)
    describe(x$beta$model)
  }

  if (!is.null(x$alpha)) {
    cli::cli_h2("Alpha diversity (linear model, marginal F tests)")
    t <- x$alpha$terms
    print(data.frame(
      metric = t$metric, term = t$term, R2 = round(t$R2, 4), F = round(t$F, 3),
      p = format.pval(t$p, digits = 2), stringsAsFactors = FALSE
    ), row.names = FALSE)
    describe(x$alpha$model)
    bad <- x$alpha$model$metric[!is.na(x$alpha$model$residual_normality_p) &
                                  x$alpha$model$residual_normality_p < 0.05]
    if (length(bad) > 0L) {
      cli::cli_alert_warning(
        "Residuals depart from normality for {paste(bad, collapse = ', ')}; treat those F-test p-values with care."
      )
    }
  }

  cli::cli_text("")
  cli::cli_alert_info(paste0(
    "Shared R2 is variance the terms explain together that no single term can claim. ",
    "Large shared R2 means the terms overlap and cannot be separated in this design; ",
    "negative shared R2 means they suppress each other."
  ))
  cli::cli_alert_info(paste0(
    "Choose these terms before fitting. If they came out of ap_screen(), say so in the methods."
  ))
  invisible(x)
}
