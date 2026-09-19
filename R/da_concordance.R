# Concordance across methods.
#
# The output people want is one list of significant features. The honest output
# is how much the methods agreed on it. A feature called by four methods and a
# feature called by one are not the same claim, and collapsing them into one
# table hides the only information that says which to trust.
#
# Agreement is computed on direction and significance. Effect magnitudes are on
# four different scales (see da.R) and cannot be compared without inventing a
# comparability that is not there.

#' Agreement between differential abundance methods
#'
#' Summarises, per feature, how many methods called it, whether they agreed on
#' direction, and which did not.
#'
#' @param da An `ap_da` result from [ap_da()].
#' @param contrast Which contrast to summarise. Defaults to the first.
#' @param min_methods Methods that must agree for a feature to enter the
#'   consensus set. Defaults to all that ran.
#'
#' @return An object of class `ap_da_concordance`: a list with `features` (one
#'   row per feature), `pairwise` (a Jaccard agreement matrix), and `consensus`.
#' @export
ap_da_concordance <- function(da, contrast = NULL, min_methods = NULL) {
  ap_assert(inherits(da, "ap_da"), "`da` must come from `ap_da()`, not {class(da)[1]}.")
  r <- da$results
  contrast <- contrast %||% r$contrast[1]
  r <- r[r$contrast == contrast, , drop = FALSE]
  ap_assert(nrow(r) > 0L, "No results for contrast `{contrast}`.")

  methods <- intersect(da$methods, unique(r$method))
  min_methods <- min_methods %||% length(methods)

  features <- sort(unique(r$feature))
  sig <- matrix(FALSE, nrow = length(features), ncol = length(methods),
                dimnames = list(features, methods))
  dir <- matrix(NA_real_, nrow = length(features), ncol = length(methods),
                dimnames = list(features, methods))
  for (m in methods) {
    sub <- r[r$method == m, ]
    idx <- match(sub$feature, features)
    sig[idx, m] <- sub$significant
    # `direction` carries the sign for structural zeros, which have no effect.
    dir[idx, m] <- if ("direction" %in% names(sub)) sub$direction else sign(sub$effect)
  }

  n_called <- rowSums(sig)
  # Direction agreement is judged only among the methods that called the
  # feature. A method that found nothing has no direction to disagree about.
  dir_called <- dir
  dir_called[!sig] <- NA_real_
  agree_dir <- apply(dir_called, 1, function(v) {
    v <- v[!is.na(v) & v != 0]
    length(v) == 0L || length(unique(v)) == 1L
  })

  feat <- data.frame(
    feature = features,
    n_methods = as.integer(n_called),
    methods = apply(sig, 1, function(v) paste(methods[v], collapse = ",")),
    direction_agrees = agree_dir,
    consensus = n_called >= min_methods & agree_dir,
    stringsAsFactors = FALSE
  )
  for (m in methods) feat[[paste0("sig_", m)]] <- sig[, m]
  feat <- feat[order(-feat$n_methods, feat$feature), ]
  rownames(feat) <- NULL

  if ("taxon_label" %in% names(r)) {
    feat$taxon_label <- r$taxon_label[match(feat$feature, r$feature)]
  }

  # Pairwise Jaccard on the significant sets.
  pairwise <- matrix(NA_real_, length(methods), length(methods),
                     dimnames = list(methods, methods))
  for (i in seq_along(methods)) {
    for (j in seq_along(methods)) {
      a <- sig[, i]; b <- sig[, j]
      union_n <- sum(a | b)
      pairwise[i, j] <- if (union_n == 0L) NA_real_ else sum(a & b) / union_n
    }
  }

  structure(
    list(features = feat, pairwise = pairwise,
         consensus = feat$feature[feat$consensus],
         n_by_method = colSums(sig),
         methods = methods, contrast = contrast,
         min_methods = min_methods, alpha = da$alpha, da = da),
    class = "ap_da_concordance"
  )
}

#' @export
print.ap_da_concordance <- function(x, ...) {
  cli::cli_h1("Method concordance: {x$contrast}")
  cli::cli_text("{length(x$methods)} methods, alpha {x$alpha}, consensus needs {x$min_methods}")

  cli::cli_h2("Called by each method")
  print(data.frame(method = names(x$n_by_method),
                   n_significant = as.integer(x$n_by_method)), row.names = FALSE)

  cli::cli_h2("Agreement (Jaccard on significant sets)")
  print(round(x$pairwise, 3))

  cli::cli_h2("How many methods called each feature")
  tab <- table(x$features$n_methods)
  print(data.frame(n_methods = names(tab), n_features = as.integer(tab)), row.names = FALSE)

  disagree <- sum(!x$features$direction_agrees)
  if (disagree > 0L) {
    cli::cli_alert_warning(
      "{disagree} feature{?s} were called by more than one method with opposite signs."
    )
  }

  cli::cli_h2("Consensus set: {length(x$consensus)} feature{?s}")
  if (length(x$consensus) > 0L) {
    top <- utils::head(x$features[x$features$consensus, ], 15L)
    cols <- intersect(c("feature", "taxon_label", "n_methods"), names(top))
    print(top[, cols], row.names = FALSE)
  }

  cli::cli_text("")
  cli::cli_alert_info(paste0(
    "The disagreement is a result. Report the consensus set and the per-method ",
    "counts together; choosing whichever method called the most is how a ",
    "differential abundance list stops being reproducible."
  ))
  invisible(x)
}
