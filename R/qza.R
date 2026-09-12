# Reading QIIME 2 artifacts without QIIME 2.
#
# A .qza is a zip archive whose single top-level directory is the artifact UUID.
# Inside it: metadata.yaml (uuid, type, format, data-size), VERSION (archive and
# framework version), checksums.sha512, a provenance/ tree, and data/ holding the
# payload in a plain format. Verified against QIIME 2 amplicon 2025.7.0,
# archive format 7.0, on 2026-09-12.
#
# Only metadata.yaml, VERSION and data/ are extracted. The provenance tree is
# large (up to 124 files in the artifacts here) and nothing downstream reads it.

#' Read the header of a QIIME 2 artifact
#'
#' Parses `metadata.yaml` and `VERSION` without extracting the payload. Useful
#' for checking what an artifact holds before committing to reading it.
#'
#' @param path Path to a `.qza` or `.qzv` file.
#'
#' @return A list with `uuid`, `type`, `format`, `archive_version`,
#'   `framework_version`, and `data_files` (the paths under `data/`).
#' @export
ap_qza_info <- function(path) {
  ap_assert(length(path) == 1L && is.character(path),
            "`path` must be a single file path.")
  ap_assert(file.exists(path), "Artifact not found: {path}")
  size <- ap_file_size(path)
  ap_assert(!is.na(size) && size > 0L, "Artifact is empty: {path}")

  listing <- tryCatch(
    utils::unzip(path, list = TRUE),
    error = function(e) {
      ap_abort("Not a readable zip archive, so not a QIIME 2 artifact: {path}\n{conditionMessage(e)}")
    }
  )
  ap_assert(nrow(listing) > 0L, "Archive is empty: {path}")

  top <- unique(sub("/.*$", "", listing$Name))
  ap_assert(
    length(top) == 1L,
    paste0("A QIIME 2 artifact has exactly one top-level directory (its UUID). ",
           "Found {length(top)}: {paste(utils::head(top, 5), collapse = ', ')}. ",
           "This is not a QIIME 2 artifact.")
  )

  meta_path <- paste0(top, "/metadata.yaml")
  ap_assert(meta_path %in% listing$Name,
            "No metadata.yaml in {path}. Not a QIIME 2 artifact.")

  con <- unz(path, meta_path)
  on.exit(close(con), add = TRUE)
  meta_lines <- readLines(con, warn = FALSE)
  meta <- ap_parse_flat_yaml(meta_lines)

  for (key in c("uuid", "type", "format")) {
    ap_assert(!is.null(meta[[key]]),
              "metadata.yaml in {path} has no `{key}` field. Artifact is malformed.")
  }
  ap_assert(
    identical(meta$uuid, top),
    "UUID mismatch in {path}: directory is `{top}` but metadata.yaml says `{meta$uuid}`."
  )

  version <- list(archive = NA_character_, framework = NA_character_)
  ver_path <- paste0(top, "/VERSION")
  if (ver_path %in% listing$Name) {
    vcon <- unz(path, ver_path)
    vlines <- readLines(vcon, warn = FALSE)
    close(vcon)
    vmeta <- ap_parse_flat_yaml(vlines)
    version$archive <- vmeta[["archive"]] %||% NA_character_
    version$framework <- vmeta[["framework"]] %||% NA_character_
  }

  data_prefix <- paste0(top, "/data/")
  data_files <- listing$Name[startsWith(listing$Name, data_prefix)]
  data_files <- data_files[!endsWith(data_files, "/")]

  list(
    uuid = meta$uuid,
    type = meta$type,
    format = meta$format,
    archive_version = version$archive,
    framework_version = version$framework,
    data_files = sub(data_prefix, "", data_files, fixed = TRUE),
    path = path
  )
}

#' Read a QIIME 2 artifact into an R object
#'
#' Reads the payload of a `.qza` and returns the natural R representation for
#' its semantic type. The QIIME 2 framework is not required; the artifact is a
#' zip archive holding data in standard formats.
#'
#' Supported semantic types:
#'
#' | Semantic type | Returns |
#' | --- | --- |
#' | `FeatureTable[*]` | integer or numeric matrix, features x samples |
#' | `Phylogeny[Rooted]`, `Phylogeny[Unrooted]` | `ape::phylo` |
#' | `FeatureData[Taxonomy]` | data frame: `feature_id`, `taxon`, `confidence` |
#' | `FeatureData[Sequence]` | `Biostrings::DNAStringSet` |
#' | `SampleData[AlphaDiversity]` | named numeric vector |
#' | `DistanceMatrix` | `stats::dist` |
#' | `PCoAResults` | list: `vectors`, `eigvals`, `prop_explained` |
#'
#' Every returned object carries the artifact's `uuid`, `type`, `format` and
#' framework version as attributes, so provenance survives the read.
#'
#' @param path Path to a `.qza` file.
#' @param tmpdir Directory to extract into. Defaults to a fresh session temp
#'   directory, removed on exit.
#'
#' @return The object described above, with attributes `qza_uuid`, `qza_type`,
#'   `qza_format`, `qza_framework`.
#' @export
ap_read_qza <- function(path, tmpdir = NULL) {
  info <- ap_qza_info(path)
  ap_assert(
    length(info$data_files) > 0L,
    "Artifact {path} has no files under data/. Nothing to read."
  )

  if (is.null(tmpdir)) {
    tmpdir <- file.path(tempdir(), paste0("amplipub-qza-", info$uuid))
    on.exit(unlink(tmpdir, recursive = TRUE, force = TRUE), add = TRUE)
  }
  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

  to_extract <- paste0(info$uuid, "/data/", info$data_files)
  utils::unzip(path, files = to_extract, exdir = tmpdir, overwrite = TRUE)
  data_dir <- file.path(tmpdir, info$uuid, "data")
  ap_assert(dir.exists(data_dir),
            "Extraction of {path} produced no data directory. Archive may be corrupt.")

  obj <- ap_parse_qza_payload(info, data_dir)

  attr(obj, "qza_uuid") <- info$uuid
  attr(obj, "qza_type") <- info$type
  attr(obj, "qza_format") <- info$format
  attr(obj, "qza_framework") <- info$framework_version
  obj
}

# Dispatch on semantic type. Kept separate from the extraction so the parsers
# can be tested against plain directories without building a zip first.
ap_parse_qza_payload <- function(info, data_dir) {
  type <- info$type
  one_file <- function(name) {
    hit <- file.path(data_dir, name)
    ap_assert(file.exists(hit),
              "Artifact of type `{type}` should contain data/{name}, and does not.")
    hit
  }

  if (startsWith(type, "FeatureTable[")) {
    return(ap_read_biom_file(one_file("feature-table.biom")))
  }
  if (startsWith(type, "Phylogeny[")) {
    return(ape::read.tree(one_file("tree.nwk")))
  }
  if (identical(type, "FeatureData[Taxonomy]")) {
    return(ap_read_q2_taxonomy(one_file("taxonomy.tsv")))
  }
  if (identical(type, "FeatureData[Sequence]")) {
    return(Biostrings::readDNAStringSet(one_file("dna-sequences.fasta")))
  }
  if (identical(type, "SampleData[AlphaDiversity]")) {
    return(ap_read_q2_alpha(one_file("alpha-diversity.tsv")))
  }
  if (identical(type, "DistanceMatrix")) {
    return(ap_read_q2_distance(one_file("distance-matrix.tsv")))
  }
  if (identical(type, "PCoAResults")) {
    return(ap_read_q2_ordination(one_file("ordination.txt")))
  }

  ap_abort(paste0(
    "AmpliPub cannot yet read QIIME 2 semantic type `{type}`.\n",
    "Files present under data/: {paste(info$data_files, collapse = ', ')}.\n",
    "Export it with `qiime tools export` and read the plain file instead."
  ))
}
