# The four backends.
#
# Each is a thin adapter: build the call the package expects, run it, and map
# its output onto the common result frame. Every inherited default that would
# change the answer is overridden explicitly here rather than accepted, and the
# reason is stated at each site.

#' @keywords internal
ap_da_formula <- function(group, covariates) {
  paste(c(group, covariates), collapse = " + ")
}

#' @keywords internal
ap_da_ancombc2 <- function(counts, meta, group, covariates, subject, alpha,
                           p_adj_method, seed) {
  ap_assert(requireNamespace("ANCOMBC", quietly = TRUE),
            "Method ancombc2 needs the ANCOMBC package.")

  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(meta, row.names = rownames(meta))
  )

  res <- ANCOMBC::ancombc2(
    data = tse,
    assay.type = "counts",
    fix_formula = ap_da_formula(group, covariates),
    rand_formula = if (is.null(subject)) NULL else paste0("(1 | ", subject, ")"),
    # `group` is required for structural-zero detection, so it is always given.
    group = group,
    # prv_cut = 0: AmpliPub filtered once already, and letting ancombc2 filter
    # again would give it a different feature set from the other methods.
    prv_cut = 0,
    lib_cut = 0,
    # ancombc2's default is "holm", far more conservative than the BH the other
    # three use. Left alone it would make ANCOM-BC2 look specific by
    # construction rather than by behaviour, so the correction is set once for
    # all four methods.
    p_adj_method = p_adj_method,
    alpha = alpha,
    # Structural zeros: a taxon absent (or nearly) from every sample of a group is
    # declared differentially abundant without a test (ANCOM-II). neg_lb adds
    # ANCOM-II's lower-bound criterion, which the help page recommends when groups
    # are large (> 30). ancombc2 drops such taxa from `res`; they are added back
    # below from `zero_ind`, or they would silently vanish from the comparison.
    struc_zero = TRUE,
    neg_lb = min(table(meta[[group]])) > 30L,
    pseudo_sens = TRUE,
    verbose = FALSE,
    n_cl = 1
  )
  out <- res$res
  ap_assert(!is.null(out) && nrow(out) > 0L, "ancombc2 returned no results.")

  lfc_cols <- grep(paste0("^lfc_", group), colnames(out), value = TRUE)
  ap_assert(length(lfc_cols) > 0L,
            "ancombc2 returned no coefficient for `{group}`.")

  zi <- res$zero_ind
  reference <- levels(meta[[group]])[1]

  do.call(rbind, lapply(lfc_cols, function(lc) {
    contrast <- sub(paste0("^lfc_", group), "", lc)
    suffix <- sub("^lfc_", "", lc)
    ss_col <- paste0("passed_ss_", suffix)
    # The sensitivity analysis re-runs ANCOM-BC2 with pseudocounts 0.1, 0.5 and 1
    # added to zeros; passed_ss says whether the call survives. ANCOMBC's own robust
    # call, diff_robust, is q < alpha AND passed_ss, and its vignette strongly
    # recommends using it for the final call. It is carried as a column so ap_da()
    # can apply it; before 2026-09-18 it was only a note, and 145 of 146 ANCOM-BC2
    # calls on Baxter were counted although they failed it.
    ap_assert(ss_col %in% colnames(out),
              "ancombc2 returned no {ss_col} column, so its calls cannot be checked for robustness.")
    passed <- as.logical(out[[ss_col]])
    note <- ifelse(passed %in% FALSE, "failed pseudocount sensitivity", NA_character_)

    estimated <- ap_da_row(
      feature = out$taxon, method = "ancombc2", contrast = contrast,
      effect = out[[lc]], effect_scale = "log fold change (natural log)",
      se = out[[paste0("se_", suffix)]],
      statistic = out[[paste0("W_", suffix)]],
      p = out[[paste0("p_", suffix)]],
      p_adj = out[[paste0("q_", suffix)]],
      note = note,
      passed_sensitivity = passed
    )
    rbind(estimated, ap_da_ancombc2_structural(zi, out$taxon, group, contrast, reference))
  }))
}

# Rows for the taxa ancombc2 set aside as structural zeros. For the contrast level vs
# reference: absent from the reference only means higher in the level (+1), absent from
# the level only means lower (-1). Absent from both, or from some other group only, gives
# no call for this contrast, and says so.
#' @keywords internal
ap_da_ancombc2_structural <- function(zi, estimated, group, level, reference) {
  if (is.null(zi) || nrow(zi) == 0L) return(NULL)
  z <- zi[!zi$taxon %in% estimated, , drop = FALSE]
  if (nrow(z) == 0L) return(NULL)
  col_of <- function(lv) paste0("structural_zero (", group, " = ", lv, ")")
  ap_assert(all(c(col_of(level), col_of(reference)) %in% names(z)),
            "ancombc2's zero_ind has no column for {level} or {reference}; cannot place its structural zeros.")
  in_lv <- as.logical(z[[col_of(level)]])
  in_ref <- as.logical(z[[col_of(reference)]])
  dirn <- ifelse(in_ref & !in_lv, 1, ifelse(in_lv & !in_ref, -1, NA_real_))
  note <- ifelse(dirn %in% 1, paste0("structural zero: absent from ", reference),
          ifelse(dirn %in% -1, paste0("structural zero: absent from ", level),
                 "structural zero in another group; not estimated for this contrast"))
  ap_da_row(
    feature = z$taxon, method = "ancombc2", contrast = level,
    effect = NA_real_, effect_scale = "log fold change (natural log)",
    se = NA_real_, statistic = NA_real_, p = NA_real_, p_adj = NA_real_,
    note = note, structural_zero = !is.na(dirn), direction = dirn
  )
}

#' @keywords internal
ap_da_aldex2 <- function(counts, meta, group, covariates, subject, reference,
                         mc_samples) {
  ap_assert(requireNamespace("ALDEx2", quietly = TRUE),
            "Method aldex2 needs the ALDEx2 package.")
  ap_assert(
    is.null(covariates) && is.null(subject),
    paste0("ALDEx2's two-group test takes no covariates or random effects. ",
           "Drop it from `method`, or run it separately without adjustment and say so.")
  )
  g <- meta[[group]]
  ap_assert(
    nlevels(g) == 2L,
    "ALDEx2 here handles two groups; `{group}` has {nlevels(g)}. Subset, or drop aldex2 from `method`."
  )

  other <- setdiff(levels(g), reference)
  res <- ALDEx2::aldex(
    reads = as.data.frame(counts),
    conditions = as.character(g),
    mc.samples = mc_samples,
    test = "t",
    effect = TRUE,
    denom = "all",
    verbose = FALSE
  )

  # ALDEx2 orders diff.btw by the sorted condition labels. Flip the sign when
  # that ordering is not the reference we were asked for, so the sign of an
  # effect means the same thing across all four methods.
  flip <- if (sort(levels(g))[1] == reference) 1 else -1

  ap_da_row(
    feature = rownames(res), method = "aldex2",
    contrast = other,
    effect = flip * res$diff.btw,
    effect_scale = "median CLR difference (log2)",
    se = res$diff.win,
    statistic = flip * res$effect,
    # Wilcoxon rather than Welch: ALDEx2's own guidance is the more conservative
    # of the two, and CLR values are not reliably normal.
    p = res$wi.ep, p_adj = res$wi.eBH,
    note = ifelse(abs(res$effect) < 1, "effect below ALDEx2's |effect| > 1 guidance", NA_character_)
  )
}

#' @keywords internal
ap_da_linda <- function(counts, meta, group, covariates, subject, alpha, p_adj_method) {
  ap_assert(requireNamespace("MicrobiomeStat", quietly = TRUE),
            "Method linda needs the MicrobiomeStat package.")

  fixed <- ap_da_formula(group, covariates)
  fml <- if (is.null(subject)) paste("~", fixed) else paste("~", fixed, "+ (1 |", subject, ")")

  res <- MicrobiomeStat::linda(
    feature.dat = counts,
    meta.dat = meta,
    formula = fml,
    feature.dat.type = "count",
    # Filtering already done once, upstream.
    prev.filter = 0,
    mean.abund.filter = 0,
    max.abund.filter = 0,
    p.adj.method = p_adj_method,
    alpha = alpha,
    verbose = FALSE
  )

  terms <- grep(paste0("^", group), names(res$output), value = TRUE)
  ap_assert(length(terms) > 0L, "linda returned no coefficient for `{group}`.")

  do.call(rbind, lapply(terms, function(tm) {
    o <- res$output[[tm]]
    ap_da_row(
      feature = rownames(o), method = "linda",
      contrast = sub(paste0("^", group), "", tm),
      effect = o$log2FoldChange, effect_scale = "log2 fold change",
      se = o$lfcSE, statistic = o$stat,
      p = o$pvalue, p_adj = o$padj
    )
  }))
}

#' @keywords internal
ap_da_maaslin2 <- function(counts, meta, group, covariates, subject, reference,
                           p_adj_method) {
  ap_assert(requireNamespace("Maaslin2", quietly = TRUE),
            "Method maaslin2 needs the Maaslin2 package.")

  # Maaslin2 writes a directory of plots and tables as a side effect. It goes to
  # a temp directory that is removed afterwards, so running a comparison does
  # not litter the working directory.
  outdir <- file.path(tempdir(), paste0("amplipub-maaslin2-", as.integer(stats::runif(1, 1e6, 9e6))))
  on.exit(unlink(outdir, recursive = TRUE, force = TRUE), add = TRUE)

  # Maaslin2 puts feature names through make.names(), which prepends "X" to any
  # name starting with a digit. ASV hashes start with a digit about 40% of the
  # time, so on real data this quietly renames a large minority of features. The
  # call still returns a full, healthy-looking result; it simply no longer joins
  # to the other methods, and the concordance is then computed across a partly
  # disjoint set. Caught on the Baxter table, where it hit 247 of 391 features.
  #
  # Renaming to safe surrogates here and mapping back afterwards makes the
  # round trip exact whatever Maaslin2 does internally.
  safe <- sprintf("apfeat%06d", seq_len(nrow(counts)))
  lookup <- stats::setNames(rownames(counts), safe)
  safe_counts <- counts
  rownames(safe_counts) <- safe

  fit <- suppressWarnings(utils::capture.output(suppressMessages(
    res <- Maaslin2::Maaslin2(
      input_data = as.data.frame(t(safe_counts)),
      input_metadata = meta,
      output = outdir,
      # Filtering already done once, upstream.
      min_abundance = 0,
      min_prevalence = 0,
      min_variance = 0,
      normalization = "TSS",
      transform = "LOG",
      analysis_method = if (is.null(subject)) "LM" else "LM",
      fixed_effects = c(group, covariates),
      random_effects = subject,
      reference = paste(group, reference, sep = ","),
      correction = p_adj_method,
      # max_significance defaults to 0.25, which is a reporting threshold, not a
      # filter; every result is kept and thresholded by AmpliPub at `alpha`.
      max_significance = 1,
      plot_heatmap = FALSE,
      plot_scatter = FALSE,
      cores = 1
    )
  )))

  out <- res$results
  ap_assert(!is.null(out) && nrow(out) > 0L, "Maaslin2 returned no results.")
  out <- out[out$metadata == group, , drop = FALSE]
  ap_assert(nrow(out) > 0L, "Maaslin2 returned no coefficient for `{group}`.")

  restored <- lookup[as.character(out$feature)]
  ap_assert(
    !anyNA(restored),
    paste0("{sum(is.na(restored))} Maaslin2 feature name{?s} could not be mapped back ",
           "to the input IDs. Maaslin2 transformed them in a way this adapter does ",
           "not reverse, and the results cannot be joined to the other methods.")
  )

  ap_da_row(
    feature = unname(restored), method = "maaslin2",
    contrast = out$value,
    effect = out$coef, effect_scale = "coefficient on log2 relative abundance",
    se = out$stderr, statistic = out$coef / out$stderr,
    p = out$pval, p_adj = out$qval
  )
}
