# Provenance.
#
# AmpliPub objects are TreeSummarizedExperiment, not a wrapper class. A wrapper
# would mean forwarding every method mia, ANCOMBC and the Bioconductor stack
# already define, and every forwarded method is a place to get it wrong. The
# provenance log rides in `S4Vectors::metadata(x)$amplipub` instead, so the
# object stays a plain TSE that any Bioconductor function accepts.
#
# Everything that transforms the data appends an entry. Seeds, thresholds,
# sample and feature counts at each boundary, and input checksums all land here,
# which is what makes a result reproducible rather than merely repeatable.

#' @keywords internal
ap_new_log <- function() {
  list(
    created = as.character(Sys.time()),
    session = ap_session_stamp(),
    inputs = list(),
    steps = list()
  )
}

#' @keywords internal
ap_session_stamp <- function() {
  pkgs <- c("AmpliPub", "vegan", "phyloseq", "mia", "TreeSummarizedExperiment",
            "ANCOMBC", "ALDEx2", "MicrobiomeStat", "Maaslin2", "ape",
            "biomformat", "ggplot2")
  versions <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)),
             error = function(e) NA_character_)
  }, character(1))
  list(
    r_version = R.version.string,
    platform = R.version$platform,
    packages = versions[!is.na(versions)]
  )
}

# md5 rather than sha512: this identifies which file was read, it is not a
# security boundary. UNC paths (\\wsl.localhost\...) warn on stat, hence the
# suppression.
#' @keywords internal
ap_checksum <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(suppressWarnings(tools::md5sum(path)))
}

#' @keywords internal
ap_log_input <- function(log, role, path, detail = list()) {
  entry <- list(
    role = role,
    path = if (is.character(path)) path else "<in-memory object>",
    bytes = if (is.character(path)) ap_file_size(path) else NA_real_,
    md5 = if (is.character(path)) ap_checksum(path) else NA_character_
  )
  log$inputs[[length(log$inputs) + 1L]] <- c(entry, detail)
  log
}

#' @keywords internal
ap_log_step <- function(x, step, detail = list(), seed = NULL) {
  log <- ap_get_log(x)
  log$steps[[length(log$steps) + 1L]] <- list(
    step = step,
    time = as.character(Sys.time()),
    seed = seed,
    detail = detail
  )
  ap_set_log(x, log)
}

#' @keywords internal
ap_get_log <- function(x) {
  if (inherits(x, "list") && !is.null(x$steps)) return(x)
  md <- S4Vectors::metadata(x)
  log <- md$amplipub
  ap_assert(
    !is.null(log),
    paste0("This object carries no AmpliPub provenance log. It was not built by ",
           "`ap_import()`. Import it through AmpliPub so the analysis is traceable.")
  )
  log
}

#' @keywords internal
ap_set_log <- function(x, log) {
  md <- S4Vectors::metadata(x)
  md$amplipub <- log
  S4Vectors::metadata(x) <- md
  x
}

#' Show what an AmpliPub object holds and how it got there
#'
#' Prints the dimensions, the metadata variables available, whether a tree and
#' taxonomy are attached, and the full provenance log: which files were read
#' with their checksums, and every step applied since, with its seed.
#'
#' @param x A `TreeSummarizedExperiment` built by [ap_import()].
#' @param steps Print the step log. Default `TRUE`.
#'
#' @return `x`, invisibly.
#' @export
ap_summary <- function(x, steps = TRUE) {
  log <- ap_get_log(x)
  counts <- SummarizedExperiment::assay(x, "counts")

  cli::cli_h1("AmpliPub object")
  cli::cli_text("{.strong {nrow(counts)}} features x {.strong {ncol(counts)}} samples")
  cli::cli_text("{.strong {format(sum(counts), big.mark = ',')}} total counts, ",
                "depth {format(min(colSums(counts)), big.mark = ',')} to ",
                "{format(max(colSums(counts)), big.mark = ',')}")

  tree <- tryCatch(TreeSummarizedExperiment::rowTree(x), error = function(e) NULL)
  cli::cli_text("Tree: {if (is.null(tree)) 'none' else paste0(length(tree$tip.label), ' tips')}")

  rd <- SummarizedExperiment::rowData(x)
  ranks <- intersect(ap_rank_names(), colnames(rd))
  cli::cli_text("Taxonomy: {if (length(ranks) == 0) 'none' else paste(ranks, collapse = ', ')}")

  cd <- SummarizedExperiment::colData(x)
  cli::cli_text("Metadata: {ncol(cd)} variable{?s}")

  cli::cli_h2("Inputs")
  for (inp in log$inputs) {
    cli::cli_li("{.field {inp$role}}: {inp$path}{if (is.na(inp$md5)) '' else paste0('  md5:', substr(inp$md5, 1, 8))}")
  }

  if (steps && length(log$steps) > 0L) {
    cli::cli_h2("Steps")
    for (st in log$steps) {
      seed_txt <- if (is.null(st$seed)) "" else paste0(" [seed ", st$seed, "]")
      cli::cli_li("{.field {st$step}}{seed_txt}: {ap_fmt_detail(st$detail)}")
    }
  }
  invisible(x)
}

#' @keywords internal
ap_fmt_detail <- function(detail) {
  if (length(detail) == 0L) return("")
  paste(
    vapply(names(detail), function(k) {
      v <- detail[[k]]
      if (length(v) > 4L) v <- c(utils::head(v, 3L), paste0("... (", length(v), " total)"))
      paste0(k, "=", paste(format(v, trim = TRUE), collapse = ","))
    }, character(1)),
    collapse = "; "
  )
}

#' Extract the provenance log
#'
#' @param x A `TreeSummarizedExperiment` built by [ap_import()].
#' @return A list with `created`, `session`, `inputs` and `steps`.
#' @export
ap_provenance <- function(x) ap_get_log(x)
