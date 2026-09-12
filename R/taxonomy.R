# Parsing taxonomy strings into ranks.
#
# Verified against Greengenes2 2024.09 output from QIIME 2 2025.7.0 on
# 2026-09-12. Two properties of real strings that a naive parser gets wrong:
#
#   1. Unassigned deep ranks are TRUNCATED, not emitted empty. A genus-level
#      call ends at `g__Blautia_A_141781` with no trailing `s__`. Splitting and
#      assuming seven fields silently shifts ranks.
#   2. GTDB names carry suffixes and numeric ids: `p__Bacillota_A_368345` is the
#      phylum the literature calls Firmicutes. AmpliPub does not rename them.
#      Renaming would hide a real disagreement between reference databases
#      behind a label that looks familiar.

#' @keywords internal
ap_rank_names <- function() {
  c("domain", "phylum", "class", "order", "family", "genus", "species")
}

#' @keywords internal
ap_rank_prefixes <- function() {
  c(domain = "d__", phylum = "p__", class = "c__", order = "o__",
    family = "f__", genus = "g__", species = "s__")
}

#' Split QIIME-style taxonomy strings into rank columns
#'
#' Accepts the semicolon-delimited, prefixed strings QIIME 2 classifiers emit
#' (`d__Bacteria; p__Bacteroidota; ...`). Ranks are assigned by their prefix,
#' never by position, so truncated strings and non-standard domain prefixes
#' (`k__` from Greengenes 13_8, `Unassigned`) land in the right column or in
#' none at all.
#'
#' Ranks not present in the string, and placeholders that carry no name
#' (`g__`, `g__uncultured`, `unidentified`), become `NA`. An empty rank and a
#' rank named "uncultured" are the same statement: not identified. Keeping them
#' apart lets an `NA`-aware aggregation quietly treat "uncultured" as a taxon
#' and pool unrelated organisms under one bar.
#'
#' @param taxa Character vector of taxonomy strings.
#' @param ranks Ranks to return. Defaults to domain through species.
#'
#' @return A data frame with one column per rank, one row per input string.
#' @export
ap_parse_taxonomy <- function(taxa, ranks = ap_rank_names()) {
  ap_assert(is.character(taxa), "`taxa` must be a character vector.")
  bad <- setdiff(ranks, ap_rank_names())
  ap_assert(
    length(bad) == 0L,
    "{cli::qty(length(bad))}Unknown rank{?s}: {paste(bad, collapse = ', ')}. Known: {paste(ap_rank_names(), collapse = ', ')}."
  )

  prefixes <- ap_rank_prefixes()
  # k__ is the Greengenes 13_8 spelling of what everything since calls d__.
  alias <- c(k__ = "d__")

  parts <- strsplit(trimws(taxa), ";", fixed = TRUE)
  out <- matrix(NA_character_, nrow = length(taxa), ncol = length(ranks),
                dimnames = list(NULL, ranks))

  for (i in seq_along(parts)) {
    fields <- trimws(parts[[i]])
    fields <- fields[nzchar(fields)]
    if (length(fields) == 0L) next
    pfx <- substr(fields, 1L, 3L)
    pfx <- ifelse(pfx %in% names(alias), alias[pfx], pfx)
    values <- trimws(substring(fields, 4L))
    keep <- pfx %in% prefixes[ranks] & nzchar(values)
    if (!any(keep)) next
    rank_of <- names(prefixes)[match(pfx[keep], prefixes)]
    out[i, rank_of] <- values[keep]
  }

  out[ap_is_unnamed_taxon(out)] <- NA_character_
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Placeholders that occupy a rank without naming anything. Matched whole, case
# insensitively, so a genuine taxon whose name merely contains "unknown" is not
# thrown away.
#' @keywords internal
ap_is_unnamed_taxon <- function(x) {
  placeholders <- c("uncultured", "unidentified", "unassigned", "unknown",
                    "uncultured bacterium", "uncultured organism",
                    "metagenome", "ambiguous_taxa", "na", "none", "-")
  v <- tolower(trimws(as.character(x)))
  res <- !is.na(v) & (v %in% placeholders | !nzchar(v))
  dim(res) <- dim(x)
  res
}

#' Build a display label for each feature at a given rank
#'
#' Returns the name at `rank`, falling back up the hierarchy when that rank is
#' unassigned, so a feature classified only to family reads
#' `f__Lachnospiraceae (unassigned genus)` rather than `NA`. Features with no
#' assignment at any rank return `"Unassigned"`.
#'
#' The fallback is marked in the label on purpose. A bar chart that silently
#' prints the family name in a genus plot is a chart that overstates its own
#' resolution.
#'
#' @param tax A data frame of rank columns, as returned by [ap_parse_taxonomy()].
#' @param rank Target rank.
#' @param mark_fallback Annotate labels that fell back. Default `TRUE`.
#'
#' @return A character vector, one label per row of `tax`.
#' @export
ap_taxon_label <- function(tax, rank = "genus", mark_fallback = TRUE) {
  ap_assert(is.data.frame(tax), "`tax` must be a data frame of rank columns.")
  ap_assert(rank %in% colnames(tax),
            "Rank `{rank}` is not a column of `tax`. Available: {paste(colnames(tax), collapse = ', ')}.")

  hierarchy <- ap_rank_names()
  target_i <- match(rank, hierarchy)
  # Ranks at or above the target, deepest first.
  candidates <- intersect(rev(hierarchy[seq_len(target_i)]), colnames(tax))
  prefixes <- ap_rank_prefixes()

  label <- rep(NA_character_, nrow(tax))
  source_rank <- rep(NA_character_, nrow(tax))
  for (r in candidates) {
    need <- is.na(label)
    if (!any(need)) break
    v <- tax[[r]][need]
    got <- !is.na(v)
    label[which(need)[got]] <- paste0(prefixes[[r]], v[got])
    source_rank[which(need)[got]] <- r
  }

  fell_back <- !is.na(source_rank) & source_rank != rank
  if (mark_fallback && any(fell_back)) {
    label[fell_back] <- paste0(label[fell_back], " (unassigned ", rank, ")")
  }
  label[is.na(label)] <- "Unassigned"
  label
}
