# Metadata assessment, run before any test is chosen.
#
# Most real studies are not the two-group comparison the default test assumes.
# They have repeated measures, a sequencing batch, and covariates that track the
# variable of interest. Deciding which test to run before looking at the
# metadata is how a paired design gets analysed as if it were independent.
#
# This stage does not choose the test. It reports what the design is, so the
# choice is made with the facts in view.

#' Assess sample metadata before choosing a test
#'
#' Classifies every metadata variable, flags likely batch or run proxies, looks
#' for repeated measures, and, when a variable of interest is named, reports
#' which other variables are associated with it and are therefore candidate
#' confounders.
#'
#' Association is screened with Kruskal-Wallis for numeric variables and a
#' chi-squared test of independence for categorical ones. These are a screen,
#' not an adjustment: a flagged variable is one to think about and report, not
#' one to regress away automatically.
#'
#' @param x A `TreeSummarizedExperiment` built by [ap_import()], or a data frame.
#' @param group Optional name of the variable of interest. When given, other
#'   variables are screened for association with it.
#' @param max_levels Above this many distinct values, a categorical variable is
#'   treated as identifier-like rather than a grouping factor. Default `20`.
#' @param force Scan even when `group` is a variable no test can use. See the
#'   section below. Default `FALSE`.
#'
#' @section The grouping variable is checked before anything downstream runs:
#' Naming a column that cannot group samples is the easiest mistake to make and
#' the hardest to see afterwards, because every stage still runs and the report
#' still renders, with confident sentences attached to a comparison that means
#' nothing. The scan therefore refuses two cases outright:
#'
#' 1. **Constant.** One distinct value, so there is no second group to compare.
#' 2. **Identifier-like.** Close to one distinct value per sample, which is a
#'    sample or subject ID rather than a factor. Every group would have one
#'    member.
#'
#' Two further cases warn but do not stop: more levels than `max_levels`, and a
#' smallest level below three samples, where within-group variance cannot be
#' estimated and the dispersion test is degenerate rather than merely weak.
#' `force = TRUE` downgrades both refusals to warnings.
#'
#' @return An object of class `ap_metadata_scan`: a list with `variables` (one
#'   row per variable), `batch_candidates`, `repeated_measures`,
#'   `confounders`, `group`, `group_status` and `group_notes`.
#' @export
ap_scan_metadata <- function(x, group = NULL, max_levels = 20L, force = FALSE) {
  meta <- ap_as_metadata_df(x)
  n <- nrow(meta)
  ap_assert(n > 0L, "Metadata has no rows.")

  if (!is.null(group)) {
    ap_assert(
      group %in% names(meta),
      "Variable `{group}` is not in the metadata. Available: {paste(names(meta), collapse = ', ')}."
    )
  }

  rows <- lapply(names(meta), function(nm) ap_describe_variable(meta[[nm]], nm, n, max_levels))
  vars <- do.call(rbind, rows)
  rownames(vars) <- NULL

  batch <- vars$variable[vars$batch_name_match & vars$role %in% c("categorical", "binary")]
  repeated <- ap_find_repeated_measures(meta, vars, max_levels)

  confounders <- NULL
  check <- list(status = NA_character_, notes = character(0))
  if (!is.null(group)) {
    check <- ap_check_group_variable(vars, group, force)
    confounders <- ap_screen_confounders(meta, group, vars)
  }

  structure(
    list(variables = vars, batch_candidates = batch,
         repeated_measures = repeated, confounders = confounders,
         group = group, group_status = check$status,
         group_notes = check$notes, n_samples = n),
    class = "ap_metadata_scan"
  )
}

# The refusal cases are the two where no test is possible at all, not the ones
# where a test is merely a bad idea. A design with 25 sites is unusual; a design
# grouped by sample ID is arithmetic nonsense. Only the second kind stops the
# run, so the guard does not block legitimate studies.
#
# Warnings are returned as well as raised. `ap_warn()` reaches the run log,
# which nobody reads afterwards; the returned notes are printed by
# `print.ap_metadata_scan()`, which is what the HTML report captures.
#' @keywords internal
ap_check_group_variable <- function(vars, group, force = FALSE) {
  r <- vars[vars$variable == group, ]
  usable <- vars$variable[vars$role %in% c("binary", "categorical") &
                            vars$variable != group]
  offer <- if (length(usable) > 0L) {
    paste0("Variables that can group these samples: ", paste(usable, collapse = ", "), ".")
  } else {
    paste0("No other variable here has between 2 and max_levels levels, so this ",
           "metadata may not describe a grouped design at all.")
  }

  fail <- function(status, headline, why) {
    msg <- c(headline,
             "i" = why,
             "i" = offer,
             "i" = "Pass `force = TRUE` to scan anyway.")
    if (!force) ap_abort(msg)
    ap_warn("`force = TRUE`: grouping variable `{group}` is {status}; downstream tests compare nothing meaningful.")
    list(status = status, notes = paste0("Forced: grouping variable is ", status, "."))
  }

  if (identical(r$role, "constant")) {
    return(fail(
      "constant",
      "Grouping variable `{group}` takes one value across every sample.",
      paste0("A comparison needs two groups. Every test downstream would be run ",
             "on a single group, and every p-value would be undefined.")
    ))
  }

  if (identical(r$role, "identifier")) {
    return(fail(
      "identifier-like",
      "Grouping variable `{group}` has {r$n_unique} distinct values, about one per sample.",
      paste0("This is a sample or subject identifier, not a factor. Each group ",
             "would hold a single sample, so any difference found would be noise. ",
             "If these are subject IDs the design has repeated measures, which ",
             "needs a mixed model rather than a group comparison.")
    ))
  }

  notes <- character(0)
  if (identical(r$role, "high-cardinality")) {
    n_lev <- r$n_unique
    ap_warn(c(
      "Grouping variable `{group}` has {n_lev} levels.",
      "i" = paste0("Tests will run, but with this many levels and few samples in ",
                   "each, the result is rarely interpretable. Check this is the ",
                   "comparison you meant.")
    ))
    notes <- c(notes, paste0(
      "Grouping variable has ", n_lev, " levels. Tests ran, but with few samples ",
      "per level the result is rarely interpretable."))
  }

  if (!is.na(r$smallest_level_n) && r$smallest_level_n < 3L) {
    small <- r$smallest_level_n
    ap_warn(c(
      "The smallest level of `{group}` holds {small} sample{?s}.",
      "i" = paste0("Below three samples the within-group variance cannot be ",
                   "estimated, so betadisper and the dispersion test are ",
                   "degenerate rather than merely underpowered. Do not report a ",
                   "dispersion p-value from this design.")
    ))
    notes <- c(notes, paste0(
      "The smallest level of the grouping variable holds ", small, " sample(s). ",
      "Within-group variance cannot be estimated below three, so the dispersion ",
      "test is degenerate and its p-value should not be reported."))
  }

  list(status = if (length(notes) == 0L) "ok" else "warned", notes = notes)
}

#' @keywords internal
ap_as_metadata_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (methods::is(x, "SummarizedExperiment")) {
    return(as.data.frame(SummarizedExperiment::colData(x)))
  }
  ap_abort("`x` must be a TreeSummarizedExperiment or a data frame, not {class(x)[1]}.")
}

# Names that usually mark a technical batch. `site` and `center` are included
# because in multi-centre studies they are the only batch proxy available when
# the sequencing run is not recorded, which is exactly the Baxter situation.
#' @keywords internal
ap_batch_name_pattern <- function() {
  "run|batch|plate|lane|flow ?cell|sequencing|seq_?run|library|extraction|site|cent(er|re)|hospital|cohort"
}

#' @keywords internal
ap_describe_variable <- function(v, name, n, max_levels) {
  n_missing <- sum(is.na(v))
  v_ok <- v[!is.na(v)]
  n_unique <- length(unique(v_ok))

  role <- if (n_unique <= 1L) {
    "constant"
  } else if (is.numeric(v)) {
    if (n_unique == 2L) "binary" else "numeric"
  } else if (n_unique == 2L) {
    "binary"
  } else if (n_unique >= n - n_missing) {
    # Every observed value distinct. That is an identifier whatever `max_levels`
    # says. Testing this against `max_levels` first, as an earlier version did,
    # meant a per-sample ID in a study with fewer samples than `max_levels` was
    # classified "categorical" and could be used as a grouping variable.
    "identifier"
  } else if (n_unique > max_levels) {
    "high-cardinality"
  } else {
    "categorical"
  }

  smallest <- if (role %in% c("categorical", "binary") && n_unique > 0L) {
    min(table(as.character(v_ok)))
  } else {
    NA_integer_
  }

  data.frame(
    variable = name,
    role = role,
    n_unique = n_unique,
    n_missing = n_missing,
    pct_missing = round(100 * n_missing / n, 1),
    smallest_level_n = smallest,
    batch_name_match = grepl(ap_batch_name_pattern(), name, ignore.case = TRUE),
    stringsAsFactors = FALSE
  )
}

# A subject identifier looks like: categorical or identifier-like, and each
# value appears more than once. That is repeated measures, and it routes the
# alpha and differential abundance tests to mixed models.
#' @keywords internal
ap_find_repeated_measures <- function(meta, vars, max_levels) {
  candidates <- character(0)
  for (nm in vars$variable) {
    v <- meta[[nm]]
    if (all(is.na(v))) next
    v_ok <- as.character(v[!is.na(v)])
    tab <- table(v_ok)
    n_unique <- length(tab)
    if (n_unique < 2L) next
    # Every value seen at least twice, and enough distinct values that this is a
    # subject id rather than an ordinary grouping factor.
    if (all(tab >= 2L) && n_unique > max_levels && max(tab) <= 20L) {
      candidates <- c(candidates, nm)
    }
  }
  candidates
}

#' @keywords internal
ap_screen_confounders <- function(meta, group, vars) {
  g <- meta[[group]]
  g_is_cat <- !is.numeric(g) || length(unique(g[!is.na(g)])) <= 10L

  testable <- vars$variable[
    vars$variable != group &
      vars$role %in% c("numeric", "binary", "categorical")
  ]

  out <- lapply(testable, function(nm) {
    v <- meta[[nm]]
    keep <- !is.na(v) & !is.na(g)
    if (sum(keep) < 10L || length(unique(v[keep])) < 2L || length(unique(g[keep])) < 2L) {
      return(NULL)
    }
    res <- tryCatch({
      if (is.numeric(v) && g_is_cat) {
        ht <- stats::kruskal.test(v[keep] ~ as.factor(as.character(g[keep])))
        list(test = "Kruskal-Wallis", statistic = unname(ht$statistic), p = ht$p.value)
      } else if (!is.numeric(v) && g_is_cat) {
        tb <- table(as.character(v[keep]), as.character(g[keep]))
        ht <- suppressWarnings(stats::chisq.test(tb))
        list(test = "Chi-squared", statistic = unname(ht$statistic), p = ht$p.value)
      } else {
        ct <- suppressWarnings(stats::cor.test(as.numeric(v[keep]), as.numeric(g[keep]),
                                               method = "spearman"))
        list(test = "Spearman", statistic = unname(ct$estimate), p = ct$p.value)
      }
    }, error = function(e) NULL)
    if (is.null(res)) return(NULL)
    data.frame(variable = nm, test = res$test, statistic = res$statistic,
               p = res$p, n = sum(keep), stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out)) return(NULL)
  out$p_adj <- stats::p.adjust(out$p, method = "BH")
  out <- out[order(out$p), ]
  rownames(out) <- NULL
  out
}

#' @export
print.ap_metadata_scan <- function(x, ...) {
  cli::cli_h1("Metadata assessment")
  cli::cli_text("{x$n_samples} samples, {nrow(x$variables)} variables")

  roles <- table(x$variables$role)
  cli::cli_text("Roles: {paste(names(roles), roles, sep = ' ', collapse = ', ')}")

  usable <- x$variables[x$variables$role %in% c("binary", "categorical"), ]
  if (nrow(usable) > 0L) {
    cli::cli_h2("Grouping variables available")
    for (i in seq_len(nrow(usable))) {
      r <- usable[i, ]
      cli::cli_li("{.field {r$variable}}: {r$n_unique} level{?s}, smallest n = {r$smallest_level_n}{if (r$n_missing > 0) paste0(', ', r$n_missing, ' missing') else ''}")
    }
  }

  if (!is.null(x$group)) {
    g <- x$variables[x$variables$variable == x$group, ]
    cli::cli_h2("Grouping variable")
    cli::cli_text(
      "{.field {x$group}}: {g$role}, {g$n_unique} level{?s}",
      "{if (!is.na(g$smallest_level_n)) paste0(', smallest n = ', g$smallest_level_n) else ''}"
    )
    if (length(x$group_notes) > 0L) {
      for (nt in x$group_notes) cli::cli_alert_warning(nt)
    }
  }

  const <- x$variables$variable[x$variables$role == "constant"]
  if (length(const) > 0L) {
    cli::cli_alert_info("Constant, cannot be tested: {paste(const, collapse = ', ')}")
  }

  high_miss <- x$variables[x$variables$pct_missing >= 20, ]
  if (nrow(high_miss) > 0L) {
    cli::cli_h2("Missingness at or above 20%")
    for (i in seq_len(nrow(high_miss))) {
      cli::cli_li("{.field {high_miss$variable[i]}}: {high_miss$pct_missing[i]}% missing")
    }
  }

  if (length(x$batch_candidates) > 0L) {
    cli::cli_h2("Batch or technical proxies")
    cli::cli_alert_warning(
      "{paste(x$batch_candidates, collapse = ', ')}"
    )
    cli::cli_text(
      "Test these alongside your variable of interest. If a technical variable ",
      "explains variation comparable to the biological one, that confound belongs ",
      "in the results, not in a limitations paragraph."
    )
  }

  if (length(x$repeated_measures) > 0L) {
    cli::cli_h2("Possible repeated measures")
    cli::cli_alert_warning("{paste(x$repeated_measures, collapse = ', ')}")
    cli::cli_text(
      "Values repeat across samples. If these are subject IDs the samples are not ",
      "independent, and tests must use a mixed model with this as a random effect."
    )
  }

  if (!is.null(x$confounders)) {
    sig <- x$confounders[x$confounders$p_adj < 0.05, ]
    cli::cli_h2("Associated with {.field {x$group}} (BH q < 0.05)")
    if (nrow(sig) == 0L) {
      cli::cli_alert_success("No variable screened as associated with {x$group}.")
    } else {
      for (i in seq_len(nrow(sig))) {
        r <- sig[i, ]
        cli::cli_li("{.field {r$variable}} ({r$test}, q = {format.pval(r$p_adj, digits = 2)})")
      }
      cli::cli_text(
        "These are candidate confounders, not established ones. A screen says the ",
        "variables move together; it does not say which adjustment is right."
      )
    }
  }
  invisible(x)
}
