# Alpha diversity.
#
# Hill numbers rather than an ad hoc index pick. Richness, Shannon and Simpson
# are not three unrelated summaries; they are one family at q = 0, 1, 2, which
# sets how strongly abundance is weighted. Reporting them together shows whether
# a difference lives in the rare taxa (q = 0) or the dominant ones (q = 2), and
# that distinction is usually the actual finding.
#
# Hill numbers are in units of effective number of species, so they are directly
# comparable across q. Shannon entropy is not: it is on a log scale whose base
# is a convention, which is why the same dataset yields 5.26 in QIIME 2 and 3.64
# in vegan. See `ap_shannon_entropy()`.

#' Alpha diversity
#'
#' Computes Hill numbers at q = 0, 1 and 2, Pielou's evenness, and Faith's
#' phylogenetic diversity when a tree is attached.
#'
#' Every value comes from an established package: richness from
#' `vegan::specnumber`, Hill q1 and q2 from `vegan::renyi(hill = TRUE)`, Shannon
#' entropy and evenness from `vegan::diversity`, Faith's PD from
#' `mia::addAlpha(index = "faith_diversity")`, Chao1 and ACE from
#' `vegan::estimateR`, and rarefaction from `vegan::rrarefy`. AmpliPub adds the
#' guards, the averaging over rarefaction iterations, and the reporting.
#'
#' @section Rarefaction:
#' Richness rises with sequencing depth, so q = 0 is not comparable across
#' samples of different depth. Rarefaction is applied by default, repeated
#' `n_iter` times and averaged, with the seed recorded in the provenance log.
#' `rarefy = FALSE` is available and is the right choice only when depths are
#' already equal, for instance when the table came from
#' `qiime diversity core-metrics` and is already rarefied.
#'
#' Rarefying is defensible here and wrong for differential abundance. AmpliPub
#' does not rarefy anywhere in the differential abundance layer.
#'
#' @section Good's coverage is not offered:
#' Good's coverage is `1 - F1/N`, built on the number of singleton features. A
#' denoiser emits no singletons by design, so on DADA2 or Deblur output the
#' statistic is exactly 1.0 for every sample regardless of how shallow the
#' sequencing was. It is not a coverage estimate on this data, it is a
#' restatement of how the table was built. [ap_goods_coverage()] computes it and
#' refuses when singletons are absent.
#'
#' @section Chao1 and ACE refuse in two situations:
#' Both estimate unseen richness from the features seen once (singletons) and
#' twice. Both are computed by `vegan::estimateR`: `chao1` is its bias-corrected
#' Chao1, `S_obs + F1(F1 - 1) / (2(F2 + 1))`, and `ace` uses its rare-abundance
#' threshold of 10. AmpliPub adds only the refusal below.
#'
#' 1. **No singletons anywhere.** Denoised tables (DADA2, Deblur) have none by
#'    design, and both estimators then equal `q0` exactly, so they would report
#'    observed richness under another name.
#' 2. **Every sample has the same depth.** That is what a rarefied table looks
#'    like. Subsampling turns counts of 2 into counts of 1 at random, so the
#'    singletons may be created by rarefaction rather than seen in the data, and
#'    the estimate grows with them.
#'
#' The check runs on the table as passed in, before any rarefaction, so
#' `rarefy = TRUE` on a denoised table still refuses. `force = TRUE` computes
#' the values anyway, with a warning that names the rule it overrode.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param metrics Which metrics to compute. Any of `"q0"`, `"q1"`, `"q2"`,
#'   `"evenness"`, `"faith_pd"`, `"shannon_entropy"`, `"chao1"`, `"ace"`.
#' @param rarefy Rarefy to a common depth before computing. Default `TRUE`.
#' @param depth Rarefaction depth. `NULL` (default) uses the smallest library
#'   size, which discards no samples. See [ap_depth_candidates()] to weigh
#'   deeper cutoffs against the samples they cost.
#' @param n_iter Rarefaction iterations to average over. Default `100`.
#' @param seed Random seed, recorded with the result.
#' @param force Compute `chao1` and `ace` even when the table has no singletons
#'   or looks rarefied. Default `FALSE`. See the section on Chao1 and ACE.
#'
#' @return An object of class `ap_alpha`: a list with `values` (a data frame of
#'   one row per sample per metric), `depth`, `n_iter`, `seed` and `dropped`.
#' @export
ap_alpha <- function(x,
                     metrics = c("q0", "q1", "q2", "evenness", "faith_pd"),
                     rarefy = TRUE,
                     depth = NULL,
                     n_iter = 100L,
                     seed = 1L,
                     force = FALSE) {

  known <- c("q0", "q1", "q2", "evenness", "faith_pd", "shannon_entropy", "chao1", "ace")
  bad <- setdiff(metrics, known)
  ap_assert(length(bad) == 0L,
            "{cli::qty(length(bad))}Unknown metric{?s}: {paste(bad, collapse = ', ')}. Known: {paste(known, collapse = ', ')}.")

  counts <- SummarizedExperiment::assay(x, "counts")
  ap_check_count_matrix(counts)

  # Before rarefaction, on purpose: the subsamples would contain singletons the
  # data never had.
  if (any(c("chao1", "ace") %in% metrics)) ap_check_richness_estimators(counts, force)

  tree <- tryCatch(TreeSummarizedExperiment::rowTree(x), error = function(e) NULL)
  if ("faith_pd" %in% metrics && is.null(tree)) {
    ap_warn(paste0(
      "Faith's PD needs a phylogeny and this object has none. Dropping it. ",
      "Pass `tree =` to `ap_import()` to get phylogenetic diversity."
    ))
    metrics <- setdiff(metrics, "faith_pd")
  }
  if (!"faith_pd" %in% metrics) tree <- NULL

  depths <- colSums(counts)
  dropped <- character(0)

  if (rarefy) {
    ap_assert(all(counts == round(counts)),
              paste0("Rarefaction subsamples reads and needs integer counts. This table ",
                     "holds fractional values. Pass `rarefy = FALSE` if it is already ",
                     "normalised, and note that richness is then not comparable across samples."))
    if (is.null(depth)) {
      depth <- min(depths)
      cli::cli_inform(paste0(
        "Rarefying to {format(depth, big.mark = ',')} reads, the smallest library. ",
        "No samples lost. `ap_depth_candidates()` shows what a deeper cutoff would cost."
      ))
    }
    ap_assert(depth >= 1, "`depth` must be at least 1, not {depth}.")
    dropped <- names(depths)[depths < depth]
    ap_assert(
      length(dropped) < ncol(counts),
      "Depth {depth} is above every library size (max {max(depths)}). Nothing would remain."
    )
    if (length(dropped) > 0L) {
      ap_warn(paste0(
        "{length(dropped)} sample{?s} below depth {depth} and dropped: ",
        "{paste(utils::head(dropped, 5), collapse = ', ')}",
        "{if (length(dropped) > 5) paste0(' and ', length(dropped) - 5, ' more') else ''}."
      ))
      counts <- counts[, setdiff(colnames(counts), dropped), drop = FALSE]
    }
    values <- ap_alpha_rarefied(counts, metrics, depth, n_iter, seed, tree)
  } else {
    if (length(unique(depths)) > 1L) {
      ap_warn(paste0(
        "`rarefy = FALSE` on a table whose depths range from {format(min(depths), big.mark = ',')} ",
        "to {format(max(depths), big.mark = ',')}. Richness (q0) rises with depth, so q0 ",
        "differences between samples may be depth differences."
      ))
    }
    n_iter <- 1L
    values <- ap_alpha_one(counts, metrics, tree)
  }

  structure(
    list(values = values,
         metrics = metrics,
         rarefied = rarefy,
         depth = if (rarefy) depth else NA_real_,
         n_iter = n_iter,
         seed = if (rarefy) seed else NA_integer_,
         dropped = dropped,
         metadata = as.data.frame(SummarizedExperiment::colData(x))),
    class = "ap_alpha"
  )
}

#' @keywords internal
ap_alpha_rarefied <- function(counts, metrics, depth, n_iter, seed, tree) {
  set.seed(seed)
  acc <- NULL
  for (i in seq_len(n_iter)) {
    sub <- ap_rarefy_matrix(counts, depth)
    one <- ap_alpha_one(sub, metrics, tree)
    acc <- if (is.null(acc)) one else {
      acc$value <- acc$value + one$value
      acc
    }
  }
  acc$value <- acc$value / n_iter
  acc
}

# Rarefaction by vegan::rrarefy: each sample drawn down to `depth` reads without
# replacement. vegan takes samples in rows, so the table is transposed in and
# back out.
#
# rrarefy warns whenever a table's smallest count is not 1, because missing
# singletons can mean the counts were multiplied. Denoised amplicon tables have
# no singletons by design, so on DADA2 or Deblur output that warning would fire
# on every one of the rarefaction iterations. It alone is silenced.
#' @keywords internal
ap_rarefy_matrix <- function(counts, depth) {
  out <- ap_muffle_warning(t(vegan::rrarefy(t(counts), depth)),
                           "function should be used for observed counts")
  dimnames(out) <- dimnames(counts)
  out
}

#' @keywords internal
ap_alpha_one <- function(counts, metrics, tree = NULL) {
  res <- list()
  sites <- t(counts)  # vegan takes samples in rows

  # Richness as an integer count. vegan::renyi gives it as exp(log S), which
  # is off by rounding error and would break exact comparisons.
  richness <- vegan::specnumber(sites)
  if ("q0" %in% metrics) res$q0 <- richness
  if (any(c("q1", "q2") %in% metrics)) {
    hill <- vegan::renyi(sites, scales = c(1, 2), hill = TRUE)
    # renyi returns a data frame for several samples and a vector for one.
    hill <- matrix(unlist(hill), nrow = nrow(sites), dimnames = list(rownames(sites), c("1", "2")))
    if ("q1" %in% metrics) res$q1 <- hill[, "1"]
    if ("q2" %in% metrics) res$q2 <- hill[, "2"]
  }
  if ("shannon_entropy" %in% metrics) res$shannon_entropy <- vegan::diversity(sites, base = 2)
  if ("evenness" %in% metrics) {
    ev <- vegan::diversity(sites) / log(richness)
    ev[richness <= 1L] <- NA_real_
    res$evenness <- ev
  }
  if ("faith_pd" %in% metrics) res$faith_pd <- ap_faith_pd_mia(counts, tree)
  if (any(c("chao1", "ace") %in% metrics)) {
    # vegan::estimateR takes sites in rows, so the table is transposed. It returns
    # one column per sample. Where ACE is undefined (every rare read a singleton)
    # vegan gives a non-finite value, which is reported as NA.
    est <- vegan::estimateR(t(counts))
    if ("chao1" %in% metrics) res$chao1 <- est["S.chao1", ]
    if ("ace" %in% metrics) {
      ace <- est["S.ACE", ]
      ace[!is.finite(ace)] <- NA_real_
      res$ace <- ace
    }
  }

  out <- do.call(rbind, lapply(names(res), function(m) {
    data.frame(sample_id = colnames(counts), metric = m,
               value = unname(res[[m]]), stringsAsFactors = FALSE)
  }))
  out[order(match(out$metric, metrics), out$sample_id), ]
}

# Chao1 and ACE extrapolate unseen richness from the features seen once. Two
# situations leave them nothing real to extrapolate from, and both are the
# normal state of amplicon data. Measured on Baxter 2016: the denoised table had
# no singletons at all, while its rarefied version had them in 448 of 487
# samples, every one created by the subsampling.
#' @keywords internal
ap_check_richness_estimators <- function(counts, force) {
  ap_assert(all(counts == round(counts)),
            "Chao1 and ACE are defined on integer counts. This table holds fractional values.")
  n <- ncol(counts)

  if (all(colSums(counts == 1) == 0)) {
    if (!force) {
      ap_abort(paste0(
        "No feature has a count of 1 in any of the {n} samples, so Chao1 and ACE would ",
        "equal observed richness (q0) exactly.\n",
        "DADA2 and Deblur remove singletons by design. The estimators then have nothing ",
        "to extrapolate from and would report q0 under another name.\n",
        "Report q0. Pass `force = TRUE` if you have a reason to want the numbers anyway."
      ))
    }
    ap_warn("`force = TRUE`: Chao1 and ACE computed on a table with no singletons, so they equal q0.")
    return(invisible("no_singletons"))
  }

  depths <- colSums(counts)
  if (n > 1L && length(unique(depths)) == 1L) {
    if (!force) {
      d <- format(depths[1], big.mark = ",")
      ap_abort(paste0(
        "Every sample has exactly {d} reads, which is what a rarefied table looks like.\n",
        "Rarefaction subsamples reads and turns counts of 2 into counts of 1 at random. ",
        "The singletons Chao1 and ACE extrapolate from may then come from the subsampling ",
        "rather than the data, and the estimate grows with them.\n",
        "Pass the unrarefied table. Pass `force = TRUE` if this table is not rarefied."
      ))
    }
    ap_warn("`force = TRUE`: Chao1 and ACE computed on an equal-depth table; its singletons may come from rarefaction.")
    return(invisible("equal_depth"))
  }
  invisible("ok")
}

#' Shannon entropy, stated in a named base
#'
#' Shannon entropy is base-dependent and the base is a convention nobody
#' states. QIIME 2 reports it in bits (base 2) because scikit-bio defaults
#' there; `vegan::diversity` reports nats (base e). The same sample gives 5.26
#' and 3.64. Neither is wrong and they are not comparable.
#'
#' Prefer the Hill number `q1`, which is `exp(H)` in nats and is base-free:
#' effective number of equally abundant species.
#'
#' Computed by `vegan::diversity(index = "shannon", base = base)`. This function
#' exists to make the base explicit and to take samples in columns.
#'
#' @param counts Numeric vector or matrix of counts, samples in columns.
#' @param base Logarithm base. Default `2`, matching QIIME 2.
#' @return Numeric vector of entropies.
#' @export
ap_shannon_entropy <- function(counts, base = 2) {
  m <- if (is.matrix(counts)) counts else matrix(counts, ncol = 1L)
  vegan::diversity(t(m), index = "shannon", base = base)
}

#' Good's coverage, with a refusal for denoised data
#'
#' Good's coverage estimates the proportion of a community's individuals that
#' belong to taxa already observed: `1 - F1/N`, where `F1` is the number of
#' features seen exactly once.
#'
#' On a denoised table this is meaningless. DADA2 and Deblur emit no singleton
#' features by construction, so `F1` is 0 and the statistic returns exactly 1.0
#' for every sample, including a sample sequenced to 500 reads. It measures how
#' the table was built, not how well the community was covered. This function
#' therefore refuses rather than returning 1.0.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param force Return the value anyway when no singletons are present.
#'   Default `FALSE`.
#' @return A named numeric vector, one value per sample.
#' @export
ap_goods_coverage <- function(x, force = FALSE) {
  counts <- SummarizedExperiment::assay(x, "counts")
  singletons <- colSums(counts == 1)
  n <- colSums(counts)
  value <- 1 - singletons / n

  if (all(singletons == 0L) && !force) {
    ap_abort(paste0(
      "No feature in this table has a count of 1 in any sample, so Good's coverage ",
      "is exactly 1.0 for all {ncol(counts)} samples.\n",
      "DADA2 and Deblur remove singletons by design, which makes the statistic a ",
      "restatement of the denoiser rather than a coverage estimate. Reporting 1.0 ",
      "would claim complete coverage for every sample including the shallowest.\n",
      "Use rarefaction curves or Hill q0 across depths instead. ",
      "Pass `force = TRUE` if you have a reason to want the number anyway."
    ))
  }
  value
}

#' What a rarefaction depth would cost
#'
#' Tabulates, across candidate depths, how many samples are retained and how
#' many reads are discarded. The depth choice is a trade between comparability
#' and sample size, and it belongs in the methods section, so it is made by
#' looking at this table rather than taken from a default.
#'
#' @param x A `TreeSummarizedExperiment` from [ap_import()].
#' @param candidates Depths to evaluate. Defaults to the deciles of the
#'   observed library sizes plus common round numbers.
#' @return A data frame: `depth`, `samples_retained`, `pct_samples_retained`,
#'   `reads_retained`, `pct_reads_retained`.
#' @export
ap_depth_candidates <- function(x, candidates = NULL) {
  depths <- colSums(SummarizedExperiment::assay(x, "counts"))
  if (is.null(candidates)) {
    q <- unique(round(stats::quantile(depths, probs = seq(0, 0.5, by = 0.05))))
    round_numbers <- c(1000, 2000, 5000, 10000, 20000, 50000)
    candidates <- sort(unique(c(q, round_numbers[round_numbers <= max(depths)])))
  }
  total_reads <- sum(depths)
  do.call(rbind, lapply(candidates, function(d) {
    keep <- depths >= d
    data.frame(
      depth = d,
      samples_retained = sum(keep),
      pct_samples_retained = round(100 * mean(keep), 1),
      reads_retained = sum(keep) * d,
      pct_reads_retained = round(100 * sum(keep) * d / total_reads, 1)
    )
  }))
}

#' @export
print.ap_alpha <- function(x, ...) {
  cli::cli_h1("Alpha diversity")
  if (x$rarefied) {
    cli::cli_text("Rarefied to {format(x$depth, big.mark = ',')} reads, ",
                  "{x$n_iter} iteration{?s} averaged, seed {x$seed}")
  } else {
    cli::cli_alert_warning("Not rarefied. Richness is not comparable across unequal depths.")
  }
  if (length(x$dropped) > 0L) {
    cli::cli_text("{length(x$dropped)} sample{?s} dropped below depth")
  }

  wide <- stats::aggregate(value ~ metric, data = x$values, FUN = function(v) {
    c(mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE),
      min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE))
  })
  summary_df <- data.frame(metric = wide$metric, round(as.data.frame(wide$value), 3))
  print(summary_df, row.names = FALSE)
  invisible(x)
}
