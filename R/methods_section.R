# A methods section, generated from what actually ran.
#
# Unlike ap_report_guide(), which describes the methods and must never vary with
# the data, this states what was done to THIS dataset: the parameters used, the
# numbers that came out, and the versions of the software that produced them.
#
# Nothing here is written from memory. Versions come from packageVersion() and
# from the conda environment export, citations come from each package's own
# CITATION file, and the numbers come from the run's own output files. A value
# that cannot be read is omitted rather than guessed.

#' A methods section for the analysis that ran
#'
#' Generates a draft methods section in past tense, naming the parameters that
#' were used, the counts that resulted, and the exact version of every tool.
#' References are taken from each package's own `CITATION` file at run time, so
#' they match the versions that actually ran and none of them are invented.
#'
#' Numbers from the upstream stages are included when `run_dir` is given and the
#' relevant summary file can be read. Anything unavailable is left out rather
#' than approximated, so the draft may need values filled in by hand.
#'
#' @param res The results list written by the workflow's analysis stage.
#' @param run_dir The run's output directory, holding `q2/`. Optional.
#' @param provenance_dir Directory holding the `env_*.yml` conda exports.
#'   Optional; without it, conda tool versions are omitted.
#'
#' @return A character vector of markdown lines.
#' @export
ap_methods_section <- function(res, run_dir = NULL, provenance_dir = NULL) {
  cfg <- res$config %||% list()
  env <- ap_env_versions(provenance_dir)
  cited <- character(0)
  cite <- function(pkg) {
    if (!pkg %in% cited) cited <<- c(cited, pkg)
    ap_cite_inline(pkg)
  }

  out <- character(0)
  add <- function(...) out <<- c(out, ..., "")
  ver <- function(name, prefix = " ") {
    v <- env[[name]]
    if (is.null(v)) "" else paste0(prefix, "v", v)
  }

  add("# Methods", "",
      "*Draft generated from the run's own parameters, outputs and software versions.",
      "Check every number against your data before submission, and add the study",
      "description, sample collection and sequencing details, which this pipeline does",
      "not know.*")

  # --- sequence processing ----------------------------------------------------
  add("## Sequence processing")

  prm <- cfg$primers
  if (!is.null(prm$forward) && !is.null(prm$reverse)) {
    add(sprintf(paste0("Primer sequences (forward %s, reverse %s) were removed with cutadapt%s ",
                       "as implemented in QIIME 2%s. Reads in which no primer was found were %s."),
                prm$forward, prm$reverse, ver("cutadapt"), ver("qiime2-version", " "),
                if (isTRUE(prm$discard_untrimmed)) "discarded" else "retained"))
  }

  q <- cfg$quality
  tl <- ap_read_run_tsv(run_dir, "q2/quality/trunc_len.tsv", header = FALSE)
  if (!is.null(q)) {
    trunc <- if (!is.null(tl)) {
      f <- ap_kv(tl, "trunc_len_f"); r <- ap_kv(tl, "trunc_len_r")
      ov <- ap_kv(tl, "expected_overlap")
      if (!is.na(f) && !is.na(r)) {
        sprintf(paste0(" Forward and reverse reads were truncated at %s and %s nt ",
                       "respectively, leaving an expected overlap of %s bp."), f, r, ov)
      } else ""
    } else ""
    add(sprintf(paste0("Truncation positions were chosen from quality profiles of %s ",
                       "subsampled reads, as the longest positions at which the median ",
                       "Phred score remained at or above %s while retaining at least %s bp ",
                       "of overlap for an expected amplicon of %s bp plus a %s bp margin.%s"),
                format(q$n %||% NA, big.mark = ","), q$min_q %||% NA,
                q$min_overlap %||% NA, q$amplicon_len %||% NA, q$margin %||% NA, trunc))
  }

  dstat <- ap_read_run_tsv(run_dir, "q2/dada2/dada2_stats.tsv")
  denoise <- sprintf(paste0("Reads were denoised, merged and chimera-filtered with DADA2%s ",
                            "as implemented in QIIME 2%s, yielding amplicon sequence variants (ASVs)."),
                     ver("bioconductor-dada2"), ver("qiime2-version", " "))
  if (!is.null(dstat) && all(c("input", "non.chimeric") %in% names(dstat))) {
    tot <- sum(dstat$input, na.rm = TRUE); keep <- sum(dstat$non.chimeric, na.rm = TRUE)
    denoise <- paste0(denoise, sprintf(
      " Of %s input read pairs, %s (%.1f%%) were retained as non-chimeric.",
      format(tot, big.mark = ","), format(keep, big.mark = ","), 100 * keep / tot))
  }
  add(denoise)

  cl <- cfg$classifier
  if (!is.null(cl$name)) {
    add(sprintf(paste0("Taxonomy was assigned with the QIIME 2 naive Bayes classifier ",
                       "(q2-feature-classifier%s) against %s, at a confidence threshold of %s."),
                ver("q2-feature-classifier", " "), cl$name, cl$confidence %||% NA))
    add(paste0("Taxon names are given as in the reference with the rank prefix removed. ",
               "Letter and numeric suffixes that separate lineages sharing a name are kept ",
               "as part of the name, and features unassigned at a rank are labelled with ",
               "the deepest rank they reached."))
  }

  add(sprintf(paste0("Representative sequences were aligned with MAFFT%s, the alignment was ",
                     "masked, and an approximately maximum-likelihood phylogeny was inferred ",
                     "with FastTree%s."), ver("mafft"), ver("fasttree")))

  cp <- cfg$collapse
  if (!is.null(cp$group_column)) {
    add(sprintf(paste0("Sequencing runs belonging to the same biological sample were pooled ",
                       "by summing counts on the `%s` column."), cp$group_column))
  }

  ft <- cfg$filter
  fs <- ap_read_run_tsv(run_dir, "q2/filtered/filter_summary.tsv")
  if (!is.null(ft)) {
    parts <- character(0)
    if (nzchar(ft$include %||% "")) {
      parts <- c(parts, sprintf("restricted to %s", ft$include))
    }
    if (nzchar(ft$exclude %||% "")) {
      ex <- trimws(strsplit(ft$exclude, ",")[[1]])
      parts <- c(parts, sprintf("with %s sequences removed", ap_and(ex)))
    }
    txt <- sprintf("Features were %s.", paste(parts, collapse = ", "))
    if (!is.null(ft$min_samples_fraction) && ft$min_samples_fraction > 0) {
      txt <- paste0(txt, sprintf(
        " Features present in fewer than %.0f%% of samples were discarded.",
        100 * ft$min_samples_fraction))
    }
    if (!is.null(fs) && nrow(fs) > 0L) {
      last <- fs[nrow(fs), ]
      txt <- paste0(txt, sprintf(
        " The filtered table held %s samples and %s features.",
        format(last$samples, big.mark = ","), format(last$features, big.mark = ",")))
    }
    add(txt)
  }

  # --- diversity --------------------------------------------------------------
  add("## Diversity")

  if (!is.null(res$depth) && !is.na(res$depth)) {
    add(sprintf(paste0("Samples were rarefied to %s reads before diversity analysis; samples ",
                       "below this depth were excluded. Rarefaction was performed under a ",
                       "fixed random seed (%s) so the subsample is reproducible."),
                format(res$depth, big.mark = ","), cfg$analysis$seed %||% NA))
  }

  if (!is.null(res$alpha_test)) {
    add(sprintf(paste0("Alpha diversity was summarised as Hill numbers of order 0, 1 and 2 ",
                       "(observed richness, the exponential of Shannon entropy, and inverse ",
                       "Simpson), Pielou's evenness, and Faith's phylogenetic diversity, ",
                       "computed with vegan (%s)."), cite("vegan")))
  }

  if (!is.null(res$permanova)) {
    mets <- res$beta$metrics %||% character(0)
    add(sprintf(paste0("Between-sample dissimilarity was measured with %s, computed with ",
                       "vegan (%s)."),
                ap_and(vapply(mets, ap_distance_label, character(1))), cite("vegan")))
  }

  # --- statistics -------------------------------------------------------------
  add("## Statistical analysis")

  an <- cfg$analysis %||% list()

  if (!is.null(res$alpha_test)) {
    tests <- unique(res$alpha_test$results$test)
    add(sprintf(paste0("Alpha diversity was compared between levels of `%s` using %s, with the ",
                       "test selected from the design and from checks of normality and equal ",
                       "variance rather than fixed in advance. Effect sizes are reported with ",
                       "95%% confidence intervals. P-values were adjusted across metrics with ",
                       "the Benjamini-Hochberg procedure."),
                an$group %||% "the grouping variable", ap_and(tests)))
  }

  if (!is.null(res$permanova)) {
    add(sprintf(paste0("Differences in composition were tested by PERMANOVA (`adonis2`, vegan; ",
                       "%s) with %s permutations, fitting each term after all others so that ",
                       "results do not depend on the order in which terms were specified. ",
                       "Homogeneity of multivariate dispersion was tested separately with ",
                       "`betadisper` and a permutation test, because PERMANOVA is sensitive to ",
                       "unequal dispersion and a significant result alone cannot distinguish a ",
                       "shift in location from a difference in spread."),
                cite("vegan"), format(an$permutations %||% 999, big.mark = ",")))
    if (!is.null(res$permanova$pairwise)) {
      add(paste0("For terms with three or more groups, pairwise PERMANOVA was run between ",
                 "each pair of groups on the distance matrix subset to that pair, with ",
                 "p-values adjusted by the Benjamini-Hochberg procedure within each term."))
    }
  }

  if (!is.null(res$ordination_diagnostics)) {
    meth <- unique(res$ordination_diagnostics$method)
    fit <- if (all(meth == "pcoa")) {
      paste0("The negative eigenvalue fraction and the variation carried by the first two ",
             "axes are reported, since a non-Euclidean distance is distorted by the projection.")
    } else if (all(meth == "nmds")) {
      "Stress is reported for each ordination as a measure of how well it represents the distances."
    } else {
      paste0("Goodness of representation is reported for each: stress for NMDS, and the ",
             "negative eigenvalue fraction with the variation on the first two axes for PCoA.")
    }
    add(sprintf("Ordinations were produced by %s (vegan; %s). %s",
                ap_and(vapply(meth, ap_ord_label, character(1))), cite("vegan"), fit))
  }

  if (!is.null(res$da)) {
    meths <- res$concordance$methods %||% unique(res$da$results$method) %||% character(0)
    labs <- vapply(meths, ap_da_label, character(1))
    keys <- vapply(meths, function(m) {
      p <- ap_da_pkg(m); if (is.na(p)) "" else cite(p)
    }, character(1))
    shown <- paste0(labs, ifelse(nzchar(keys), paste0(" (", keys, ")"), ""))
    primary <- an$da_primary %||% "ancombc2"
    add(sprintf(paste0("Differential abundance was assessed with %s, all run on the same ",
                       "pre-filtered feature set%s.%s%s Because these methods differ in their ",
                       "assumptions about compositionality, zeros and normalization, the other ",
                       "methods are reported alongside it, and a consensus set was defined as ",
                       "features called by at least %s methods that also agreed on the direction ",
                       "of the effect."),
                ap_and(shown),
                if (!is.null(an$da_reference)) sprintf(" against the reference level `%s`",
                                                       an$da_reference) else "",
                if ("ancombc2" %in% meths) paste0(
                  " ANCOM-BC2 calls were required to pass its pseudocount sensitivity ",
                  "analysis (the package's robust call).") else "",
                if (primary %in% meths) sprintf(paste0(
                  " %s was designated the primary method in the analysis configuration; its ",
                  "effect estimates are reported with 95%% Wald intervals (estimate %s 1.96 ",
                  "standard errors)."), ap_da_method_title(primary), "\u00b1") else "",
                res$concordance$min_methods %||% length(meths)))
  }

  if (!is.null(res$screen)) {
    add(sprintf(paste0("An exploratory screen tested every usable metadata variable against ",
                       "alpha and beta diversity. Results were ranked by variance explained ",
                       "adjusted for degrees of freedom rather than by p-value, because ",
                       "unadjusted R2 increases with the number of levels a variable has. A ",
                       "single Benjamini-Hochberg correction was applied across all tests, and ",
                       "rank stability was estimated from %s resamples. These results are ",
                       "hypothesis-generating and are not reported as findings."),
                format(an$n_resample %||% 100, big.mark = ",")))
  }

  if (!is.null(res$explains)) {
    add(paste0("Terms that the screen identified were then fitted together in a single model, ",
               "with each term's contribution assessed after all others (marginal sums of ",
               "squares). Variance that could not be attributed to any single term is reported ",
               "separately."))
  }

  if (!is.null(res$normalization)) {
    nrm <- unlist(an$normalizations %||% character(0))
    labs <- vapply(nrm, ap_norm_label, character(1))
    if ("css" %in% nrm) {
      labs[nrm == "css"] <- sprintf("cumulative-sum scaling (%s)", cite("metagenomeSeq"))
    }
    add(sprintf(paste0("To establish whether conclusions depended on the normalization chosen, ",
                       "the comparison was repeated under %s. Agreement across normalizations ",
                       "is reported."), ap_and(labs)))
  }

  # --- software ---------------------------------------------------------------
  add("## Software and reproducibility")

  sw <- sprintf("Analyses were performed in R %s using AmpliPub %s",
                paste(R.version$major, R.version$minor, sep = "."),
                as.character(utils::packageVersion("AmpliPub")))
  if (!is.null(env[["qiime2-version"]])) {
    sw <- paste0(sw, sprintf(", with upstream processing in QIIME 2 %s", env[["qiime2-version"]]))
  }
  sw <- paste0(sw, sprintf(paste0(". All random processes used a fixed seed (%s). Exact package ",
                                  "versions, the software environment and the analysis code ",
                                  "revision are recorded in the provenance directory accompanying ",
                                  "this report."), an$seed %||% NA))
  add(sw)

  # --- references -------------------------------------------------------------
  refs <- ap_reference_list(cited)
  if (length(refs)) add(c("## References", "", refs))

  out
}

#' @keywords internal
ap_env_versions <- function(provenance_dir) {
  if (is.null(provenance_dir) || !dir.exists(provenance_dir)) return(list())
  ymls <- list.files(provenance_dir, pattern = "^env_.*\\.yml$", full.names = TRUE)
  if (!length(ymls)) return(list())
  lines <- unlist(lapply(ymls, function(f) tryCatch(readLines(f, warn = FALSE),
                                                    error = function(e) character(0))))
  dep <- grep("^\\s*-\\s*[A-Za-z0-9_.-]+=", lines, value = TRUE)
  nm <- sub("^\\s*-\\s*([A-Za-z0-9_.-]+)=.*$", "\\1", dep)
  vr <- sub("^\\s*-\\s*[A-Za-z0-9_.-]+=([^=]+).*$", "\\1", dep)
  keep <- !duplicated(nm)
  vs <- stats::setNames(as.list(vr[keep]), nm[keep])
  # QIIME 2 itself is versioned by its plugins, which all share the release tag.
  q2 <- vs[["q2-diversity"]] %||% vs[["q2-dada2"]] %||% vs[["q2-types"]]
  if (!is.null(q2)) vs[["qiime2-version"]] <- sub("\\.0$", "", q2)
  vs
}

# `header` matters: the key/value summaries such as trunc_len.tsv have no header
# row, and reading one with header = TRUE silently consumes the first record.
# That is how the truncation lengths went missing from the first draft of this.
#' @keywords internal
ap_read_run_tsv <- function(run_dir, rel, header = TRUE) {
  if (is.null(run_dir)) return(NULL)
  p <- file.path(run_dir, rel)
  if (!file.exists(p)) return(NULL)
  tryCatch(utils::read.delim(p, comment.char = "#", header = header,
                             stringsAsFactors = FALSE),
           error = function(e) NULL)
}

# Several upstream summaries are two-column key/value tables.
#' @keywords internal
ap_kv <- function(df, key) {
  if (is.null(df) || ncol(df) < 2L) return(NA_character_)
  hit <- which(df[[1]] == key)
  if (!length(hit)) return(NA_character_)
  as.character(df[[2]][hit[1]])
}

#' @keywords internal
ap_and <- function(x) {
  x <- x[nzchar(x)]
  if (!length(x)) return("")
  if (length(x) == 1L) return(x)
  paste0(paste(utils::head(x, -1L), collapse = ", "), " and ", utils::tail(x, 1L))
}

#' @keywords internal
ap_da_label <- function(m) {
  switch(m, ancombc2 = "ANCOM-BC2", aldex2 = "ALDEx2", linda = "LinDA",
         maaslin2 = "MaAsLin2", m)
}

#' @keywords internal
ap_da_pkg <- function(m) {
  switch(m, ancombc2 = "ANCOMBC", aldex2 = "ALDEx2", linda = "MicrobiomeStat",
         maaslin2 = "Maaslin2", NA_character_)
}

#' @keywords internal
ap_distance_label <- function(m) {
  switch(m,
         bray_curtis = "Bray-Curtis dissimilarity",
         jaccard = "Jaccard distance",
         unweighted_unifrac = "unweighted UniFrac",
         weighted_unifrac = "weighted UniFrac",
         aitchison = "Aitchison distance",
         gsub("_", " ", m))
}

#' @keywords internal
ap_ord_label <- function(m) {
  switch(m,
         pcoa = "principal coordinates analysis (PCoA)",
         nmds = "non-metric multidimensional scaling (NMDS)",
         m)
}

#' @keywords internal
ap_norm_label <- function(n) {
  switch(n,
         tss = "total-sum scaling",
         css = "cumulative-sum scaling",
         clr = "a centred log-ratio transform with Aitchison distance",
         rarefy = "rarefaction",
         n)
}

# Author and year are read from the package's own CITATION. When that cannot be
# parsed the package name and version are used instead, because inventing a
# reference is worse than printing a less elegant one.
#' @keywords internal
ap_cite_inline <- function(pkg) {
  fallback <- function() {
    v <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NULL)
    if (is.null(v)) pkg else paste0(pkg, " ", v)
  }
  tryCatch({
    b <- utils::citation(pkg)[[1]]
    au <- b$author
    yr <- b$year
    if (is.null(au) || !length(au) || is.null(yr)) return(fallback())
    fam <- format(au[[1]], include = "family")
    if (!nzchar(fam)) return(fallback())
    paste0(fam, if (length(au) > 1L) " et al." else "", ", ", yr)
  }, error = function(e) fallback())
}

# Some CITATION files contain literal asterisks and encoded spaces that survive
# formatting and render as noise. Those are stripped. Nothing else about the
# entry is altered: authors, year, DOI and URL are exactly what the package
# supplied, because a reference nobody can verify is worse than an ugly one.
#' @keywords internal
ap_clean_reference <- function(txt) {
  s <- paste(txt, collapse = " ")
  # Remove each artefact together with the whitespace it introduced, rather than
  # stripping the asterisks and then trying to repair spacing afterwards. A
  # general "close the gap after a quote" rule also closes the legitimate gap
  # after a closing quote, which runs the title into the journal name.
  # Whitespace is collapsed first: the artefact is followed by a line break as
  # often as by a space, and removing it before the break becomes a space leaves
  # the gap behind.
  s <- gsub("\\s+", " ", s)
  s <- gsub("\\*{2,}%20", "", s)
  s <- gsub("\\*{2,} *", "", s)
  s <- gsub("%20", " ", s, fixed = TRUE)
  trimws(s)
}

#' @keywords internal
ap_reference_list <- function(pkgs) {
  if (!length(pkgs)) return(character(0))
  entries <- vapply(unique(pkgs), function(p) {
    tryCatch({
      paste0("- ", ap_clean_reference(format(utils::citation(p)[[1]], style = "text")))
    }, error = function(e) {
      v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e2) "")
      paste0("- ", p, if (nzchar(v)) paste0(" version ", v) else "",
             ". Citation unavailable from the installed package; add it by hand.")
    })
  }, character(1))
  entries <- sort(entries)
  c(entries, "",
    "*References were generated from the CITATION file of each installed package, so",
    "they match the versions that ran. Verify them against the journal's style before",
    "submission.*")
}
