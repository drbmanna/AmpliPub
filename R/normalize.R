# Normalization, and whether a conclusion survives changing it.
#
# There is no agreed normalization for amplicon counts. Rarefying is defensible for
# diversity and wrong for differential abundance; total-sum scaling ignores
# compositionality; CLR needs a zero replacement; CSS rests on a quantile choice. Picking
# one silently makes the result depend on a choice nobody reported. So the conclusion is
# computed under each, and the table says whether it held.
#
# Each normalization is a package call. AmpliPub adds the zero-only replacement before the
# CLR, which vegan's own pseudocount does not do, and the comparison around them.

#' Normalize a count table
#'
#' | Method | Computed by |
#' | --- | --- |
#' | `tss` | `vegan::decostand(method = "total")`: proportions |
#' | `css` | `metagenomeSeq::cumNorm` with `cumNormStatFast`, then `MRcounts(norm = TRUE)` |
#' | `clr` | zeros replaced by `pseudocount`, then `vegan::decostand(method = "clr")` |
#' | `rarefy` | `vegan::rrarefy` to `depth`, samples below it dropped |
#'
#' @section The CLR zero replacement:
#' `vegan::decostand(method = "clr", pseudocount = )` adds the pseudocount to every
#' count, which makes the result change when a sample's depth is scaled. Here only zeros
#' are replaced, then the CLR is taken with no pseudocount, so a sample with no zeros is
#' unaffected by its depth.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param method One of `"tss"`, `"css"`, `"clr"`, `"rarefy"`.
#' @param depth Rarefaction depth for `"rarefy"`. `NULL` uses the smallest library.
#' @param seed Random seed for `"rarefy"`.
#' @param pseudocount Value substituted for zeros before the CLR. Default `0.5`.
#' @param css_p Quantile for CSS. `NULL` lets `metagenomeSeq::cumNormStatFast` choose it.
#' @param css_scale The `sl` scaling passed to `metagenomeSeq::MRcounts`. Default `1000`,
#'   metagenomeSeq's own default, passed explicitly.
#'
#' @return A features-by-samples numeric matrix, with the settings used in its attributes
#'   (`method`, and `css_p`, `depth` or `pseudocount` where they apply).
#' @export
ap_normalize <- function(x,
                         method = c("tss", "css", "clr", "rarefy"),
                         depth = NULL,
                         seed = 1L,
                         pseudocount = 0.5,
                         css_p = NULL,
                         css_scale = 1000) {
  method <- match.arg(method)
  counts <- SummarizedExperiment::assay(x, "counts")
  ap_check_count_matrix(counts)

  out <- switch(
    method,
    tss = t(vegan::decostand(t(counts), method = "total")),
    clr = {
      ap_assert(pseudocount > 0, "`pseudocount` must be positive; the CLR is undefined at zero.")
      m <- counts
      m[m == 0] <- pseudocount
      res <- t(vegan::decostand(t(m), method = "clr"))
      attr(res, "pseudocount") <- pseudocount
      res
    },
    css = {
      ap_assert(
        requireNamespace("metagenomeSeq", quietly = TRUE),
        "CSS normalization needs the metagenomeSeq package. Install it, or leave `css` out of the comparison."
      )
      obj <- metagenomeSeq::newMRexperiment(counts)
      p <- css_p %||% metagenomeSeq::cumNormStatFast(obj)
      obj <- metagenomeSeq::cumNorm(obj, p = p)
      res <- metagenomeSeq::MRcounts(obj, norm = TRUE, log = FALSE, sl = css_scale)
      attr(res, "css_p") <- p
      attr(res, "css_scale") <- css_scale
      res
    },
    rarefy = {
      ap_assert(all(counts == round(counts)), "Rarefaction needs integer counts.")
      depth <- depth %||% min(colSums(counts))
      keep <- colSums(counts) >= depth
      if (any(!keep)) {
        ap_warn("{sum(!keep)} sample{?s} below depth {depth} dropped before rarefying.")
      }
      set.seed(seed)
      res <- ap_rarefy_matrix(counts[, keep, drop = FALSE], depth)
      attr(res, "depth") <- depth
      attr(res, "seed") <- seed
      res
    }
  )
  attr(out, "method") <- method
  out
}

#' Does a conclusion survive a change of normalization?
#'
#' Normalizes the table each way, computes a distance suited to that normalization,
#' runs [ap_permanova()] with its dispersion check for each term, and reports whether the
#' verdict for each term is the same under every normalization.
#'
#' | Normalization | Distance |
#' | --- | --- |
#' | `tss`, `css`, `rarefy` | Bray-Curtis, `vegan::vegdist(method = "bray")` |
#' | `clr` | Euclidean on the CLR, which is Aitchison distance |
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param terms Metadata variables to test, as for [ap_permanova()].
#' @param methods Normalizations to compare. Default all four.
#' @param depth Rarefaction depth for `"rarefy"`.
#' @param pseudocount Zero replacement for `"clr"`.
#' @param permutations Permutations for each PERMANOVA. Default `999`.
#' @param seed Random seed, recorded with the result.
#'
#' @return An object of class `ap_normalization_sensitivity`: `results` (one row per
#'   normalization per term) and `agreement` (one row per term, with the verdicts and
#'   whether they are all the same).
#' @export
ap_normalization_sensitivity <- function(x,
                                         terms,
                                         methods = c("tss", "css", "clr", "rarefy"),
                                         depth = NULL,
                                         pseudocount = 0.5,
                                         permutations = 999L,
                                         seed = 1L) {
  methods <- unique(match.arg(methods, c("tss", "css", "clr", "rarefy"), several.ok = TRUE))
  ap_assert(length(methods) >= 2L,
            "A sensitivity analysis needs at least two normalizations to compare.")
  if (inherits(terms, "formula")) terms <- all.vars(terms)
  meta <- as.data.frame(SummarizedExperiment::colData(x))

  rows <- list()
  for (mth in methods) {
    m <- ap_normalize(x, method = mth, depth = depth, seed = seed, pseudocount = pseudocount)
    distance <- if (mth == "clr") "aitchison" else "bray_curtis"
    d <- if (mth == "clr") stats::dist(t(m)) else vegan::vegdist(t(m), method = "bray")
    beta <- structure(
      list(distances = stats::setNames(list(d), distance), metrics = distance,
           rarefied = mth == "rarefy", depth = attr(m, "depth") %||% NA_real_,
           seed = seed, dropped = character(0), pseudocount = pseudocount,
           n_samples = ncol(m), metadata = meta[colnames(m), , drop = FALSE]),
      class = "ap_beta"
    )
    pn <- suppressMessages(ap_permanova(beta, terms, permutations = permutations, seed = seed))
    res <- merge(pn$results, pn$dispersion[, c("metric", "term", "dispersion_p")],
                 by = c("metric", "term"), all.x = TRUE)
    res <- merge(res, pn$interpretation[, c("metric", "term", "verdict")],
                 by = c("metric", "term"), all.x = TRUE)
    rows[[mth]] <- data.frame(normalization = mth, distance = distance, term = res$term,
                              n = res$n, R2 = res$R2, pseudo_F = res$pseudo_F, p = res$p,
                              dispersion_p = res$dispersion_p, verdict = res$verdict,
                              stringsAsFactors = FALSE)
  }
  results <- do.call(rbind, rows)
  rownames(results) <- NULL

  agreement <- do.call(rbind, lapply(terms, function(tm) {
    r <- results[results$term == tm, ]
    data.frame(term = tm,
               verdicts = paste(sprintf("%s: %s", r$normalization, r$verdict), collapse = "; "),
               same_verdict = length(unique(r$verdict)) == 1L,
               significant_under_all = all(r$p < 0.05),
               significant_under_any = any(r$p < 0.05),
               stringsAsFactors = FALSE)
  }))

  structure(list(results = results, agreement = agreement, methods = methods, terms = terms,
                 permutations = permutations, seed = seed, pseudocount = pseudocount),
            class = "ap_normalization_sensitivity")
}

#' @export
print.ap_normalization_sensitivity <- function(x, ...) {
  cli::cli_h1("Normalization sensitivity")
  cli::cli_text("{length(x$methods)} normalizations: {paste(x$methods, collapse = ', ')} | ",
                "{x$permutations} permutations | seed {x$seed}")
  r <- x$results
  print(data.frame(normalization = r$normalization, distance = r$distance, term = r$term,
                   R2 = round(r$R2, 4), p = format.pval(r$p, digits = 2),
                   disp_p = format.pval(r$dispersion_p, digits = 2), verdict = r$verdict,
                   stringsAsFactors = FALSE), row.names = FALSE)
  cli::cli_h2("Does the conclusion hold?")
  for (i in seq_len(nrow(x$agreement))) {
    a <- x$agreement[i, ]
    if (a$same_verdict) {
      cli::cli_alert_success("{.field {a$term}}: the same verdict under every normalization.")
    } else {
      cli::cli_alert_warning(
        "{.field {a$term}}: the verdict depends on the normalization. Report it that way."
      )
    }
  }
  invisible(x)
}
