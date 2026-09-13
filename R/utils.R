#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

# Every guard in AmpliPub routes through these two, so failure messages have one
# shape and one place to change. `message` is interpolated by cli, so `{x}` in a
# message picks up `x` from the calling frame.

#' @keywords internal
ap_abort <- function(message, ..., .envir = parent.frame()) {
  cli::cli_abort(message, ..., .envir = .envir, call = NULL)
}

#' @keywords internal
ap_assert <- function(cond, message, ..., .envir = parent.frame()) {
  if (!isTRUE(cond)) ap_abort(message, ..., .envir = .envir)
  invisible(TRUE)
}

#' @keywords internal
ap_warn <- function(message, ..., .envir = parent.frame()) {
  cli::cli_warn(message, ..., .envir = .envir)
}

# Some package calls emit a warning that is known and expected on amplicon data
# every time they run. A warning that always fires teaches people to ignore
# warnings, so exactly that message is silenced and every other one still shows.
#' @keywords internal
ap_muffle_warning <- function(expr, pattern) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(pattern, conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}

# metadata.yaml and VERSION inside a QIIME 2 artifact are flat `key: value`
# documents, four keys at most, no nesting and no lists. A full YAML parser is a
# dependency we do not need for that. Anything nested is rejected rather than
# silently half-parsed.
#' @keywords internal
ap_parse_flat_yaml <- function(lines) {
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[!startsWith(trimws(lines), "#")]
  if (length(lines) == 0L) return(list())

  indented <- grepl("^[[:space:]]", lines)
  if (any(indented)) {
    ap_abort(paste0(
      "Expected a flat `key: value` document and found indentation, ",
      "which means nesting this parser does not handle."
    ))
  }

  out <- list()
  for (ln in lines) {
    m <- regexec("^([^:]+):[[:space:]]*(.*)$", ln)
    parts <- regmatches(ln, m)[[1]]
    if (length(parts) != 3L) next
    out[[trimws(parts[2])]] <- trimws(parts[3])
  }
  out
}

# BIOM comes in two incompatible encodings: v1.0 is JSON, v2.x is HDF5.
# `biomformat::read_biom` tries JSON first and falls back, which means reading a
# v2.1 file emits a JSON lexer warning about the HDF5 magic bytes on every call.
# A warning that always fires is a warning nobody reads, so the encoding is
# detected from the magic bytes and the right reader is called directly.
#' @keywords internal
ap_read_biom_file <- function(file) {
  magic <- readBin(file, what = "raw", n = 4L)
  hdf5_magic <- as.raw(c(0x89, 0x48, 0x44, 0x46))  # \x89 H D F
  b <- if (length(magic) == 4L && identical(magic, hdf5_magic)) {
    biomformat::biom(biomformat::read_hdf5_biom(file))
  } else {
    biomformat::read_biom(file)
  }
  m <- methods::as(biomformat::biom_data(b), "matrix")
  storage.mode(m) <- "double"
  m
}

# file.info() on a UNC path (\\wsl.localhost\...) warns that it cannot resolve
# the file owner on every call. The size is all we want.
#' @keywords internal
ap_file_size <- function(path) {
  suppressWarnings(file.size(path))
}

# QIIME 2 taxonomy.tsv: header `Feature ID`, `Taxon`, and `Confidence` when the
# classifier emitted one. Renamed to syntactic names so downstream code never
# needs backticks.
#' @keywords internal
ap_read_q2_taxonomy <- function(file) {
  df <- utils::read.delim(file, sep = "\t", header = TRUE, quote = "",
                          check.names = FALSE, stringsAsFactors = FALSE,
                          colClasses = "character")
  ap_assert(ncol(df) >= 2L,
            "Taxonomy file {file} has {ncol(df)} column(s); expected at least Feature ID and Taxon.")
  names(df)[1:2] <- c("feature_id", "taxon")
  if (ncol(df) >= 3L && grepl("^confidence$", names(df)[3], ignore.case = TRUE)) {
    names(df)[3] <- "confidence"
    df$confidence <- suppressWarnings(as.numeric(df$confidence))
  }
  ap_assert(!anyDuplicated(df$feature_id),
            "Taxonomy file {file} has duplicate feature IDs.")
  df
}

#' @keywords internal
ap_read_q2_alpha <- function(file) {
  df <- utils::read.delim(file, sep = "\t", header = TRUE, quote = "",
                          check.names = FALSE, stringsAsFactors = FALSE)
  ap_assert(ncol(df) == 2L,
            "Alpha diversity file {file} has {ncol(df)} columns; expected exactly 2.")
  out <- suppressWarnings(as.numeric(df[[2]]))
  names(out) <- as.character(df[[1]])
  attr(out, "metric") <- names(df)[2]
  out
}

#' @keywords internal
ap_read_q2_distance <- function(file) {
  df <- utils::read.delim(file, sep = "\t", header = TRUE, row.names = 1,
                          quote = "", check.names = FALSE)
  m <- as.matrix(df)
  ap_assert(nrow(m) == ncol(m),
            "Distance matrix {file} is {nrow(m)} x {ncol(m)}; a distance matrix must be square.")
  ap_assert(identical(rownames(m), colnames(m)),
            "Distance matrix {file} has row and column IDs in different order or content.")
  diag_max <- max(abs(diag(m)))
  ap_assert(diag_max < 1e-8,
            "Distance matrix {file} has a non-zero diagonal (max |d(i,i)| = {signif(diag_max, 3)}).")
  stats::as.dist(m)
}

# scikit-bio OrdinationResults text format. Sections are separated by a blank
# line; each opens with `Name<TAB>nrow[<TAB>ncol]`. Verified against QIIME 2
# 2025.7.0 output on 2026-09-12.
#' @keywords internal
ap_read_q2_ordination <- function(file) {
  lines <- readLines(file, warn = FALSE)
  header_idx <- grep("^[A-Za-z][A-Za-z ]*\t", lines)
  ap_assert(length(header_idx) > 0L,
            "No section headers found in ordination file {file}.")

  section <- function(name) {
    hit <- header_idx[startsWith(lines[header_idx], paste0(name, "\t"))]
    if (length(hit) == 0L) return(NULL)
    hdr <- strsplit(lines[hit[1]], "\t", fixed = TRUE)[[1]]
    n <- suppressWarnings(as.integer(hdr[2]))
    if (is.na(n) || n == 0L) return(NULL)
    body <- lines[(hit[1] + 1L):(hit[1] + n)]
    strsplit(body, "\t", fixed = TRUE)
  }

  eig_raw <- section("Eigvals")
  ap_assert(!is.null(eig_raw), "Ordination file {file} has no eigenvalues.")
  eigvals <- as.numeric(eig_raw[[1]])

  prop_raw <- section("Proportion explained")
  prop <- if (is.null(prop_raw)) eigvals / sum(eigvals) else as.numeric(prop_raw[[1]])

  site_raw <- section("Site")
  ap_assert(!is.null(site_raw), "Ordination file {file} has no sample coordinates.")
  ids <- vapply(site_raw, `[`, character(1), 1L)
  coords <- t(vapply(site_raw, function(r) as.numeric(r[-1]), numeric(length(site_raw[[1]]) - 1L)))
  rownames(coords) <- ids
  colnames(coords) <- paste0("PC", seq_len(ncol(coords)))

  ap_assert(
    ncol(coords) == length(eigvals),
    paste0("Ordination file {file} is internally inconsistent: {ncol(coords)} coordinate ",
           "axes but {length(eigvals)} eigenvalues.")
  )

  list(vectors = coords, eigvals = eigvals, prop_explained = prop)
}
