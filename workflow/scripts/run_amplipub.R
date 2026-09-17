# AmpliPub analysis stage of the workflow. Called by Snakemake's `script:` directive,
# which provides the `snakemake` object (inputs, outputs, config, threads, log).
#
# What runs on which table, and why:
# - Alpha and beta diversity, PERMANOVA, the screen and ap_explains use QIIME 2's rarefied
#   table from 10_diversity. That is the table the recorded Baxter results came from, so
#   those numbers can be checked against the earlier runs.
# - Composition, differential abundance and the normalization comparison use the filtered,
#   unrarefied table. Differential abundance must not be rarefied, and the normalization
#   comparison needs the raw counts to normalize.
#
# Every table goes to tables/, every figure to figures/, and everything needed to trace the
# run to provenance/.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con)
sink(log_con, type = "message")

suppressPackageStartupMessages(library(AmpliPub))

inp <- snakemake@input
cfg <- snakemake@config
an <- cfg$analysis
out_dir <- dirname(snakemake@output[["results"]])
dir_tables <- file.path(out_dir, "tables")
dir_figures <- file.path(out_dir, "figures")
dir_prov <- file.path(out_dir, "provenance")
for (d in c(dir_tables, dir_figures, dir_prov)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

started <- Sys.time()
set.seed(an$seed)

write_tsv <- function(df, name) {
  utils::write.table(df, file.path(dir_tables, paste0(name, ".tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
save_plot <- function(p, name, width = 7, height = 5) {
  ggplot2::ggsave(file.path(dir_figures, paste0(name, ".png")), p,
                  width = width, height = height, dpi = 300)
}
# A figure that cannot be drawn (for example, no feature called by any DA method) is a
# result, not a crash. It is logged and the run continues.
try_plot <- function(expr, name, ...) {
  p <- tryCatch(expr, error = function(e) {
    message("figure ", name, " not drawn: ", conditionMessage(e))
    NULL
  })
  if (!is.null(p)) save_plot(p, name, ...)
  invisible(p)
}
section <- function(title) cat("\n==================== ", title, " ====================\n", sep = "")

# --- inputs ---------------------------------------------------------------------------------

section("inputs")
read_metadata <- function(path) {
  m <- utils::read.delim(path, na.strings = c("", "NA"), check.names = FALSE,
                         stringsAsFactors = FALSE)
  rownames(m) <- as.character(m[["sample-id"]])
  m[["sample-id"]] <- NULL
  m
}
meta <- read_metadata(inp[["metadata"]])
tree <- ap_read_qza(inp[["tree"]])
taxonomy <- ap_read_qza(inp[["taxonomy"]])

x <- ap_import(ap_read_qza(inp[["table"]]), meta, tree = tree, taxonomy = taxonomy)
ap_summary(x)

settings <- utils::read.delim(inp[["diversity_settings"]], stringsAsFactors = FALSE)
depth <- as.integer(settings$value[settings$setting == "sampling_depth"])
stopifnot(length(depth) == 1L, !is.na(depth))
cat("rarefaction depth from 10_diversity:", depth, "\n")

# Rarefy in R, under this run's seed. Decided 2026-09-16, after the same config run twice
# gave different alpha, beta, PERMANOVA and screen numbers: QIIME 2's
# `core-metrics-phylogenetic` subsamples once from an RNG it does not expose. It has no
# `--p-random-seed` (unlike `feature-table rarefy`), and its help says the seed "Defaults to
# a random seed". Recomputing here also follows the locked design decision that everything
# is recomputed in R and QIIME 2 is the cross-check, never the source.
rarefied_counts <- ap_normalize(x, method = "rarefy", depth = depth, seed = an$seed)
x_rarefied <- x[rownames(rarefied_counts), colnames(rarefied_counts)]
SummarizedExperiment::assay(x_rarefied, "counts") <- rarefied_counts
stopifnot(length(unique(colSums(SummarizedExperiment::assay(x_rarefied, "counts")))) == 1L)
cat("rarefied in R:", ncol(x_rarefied), "samples at depth", depth,
    "with seed", an$seed, "\n")

# QIIME 2's own rarefied table, kept as the independent cross-check it is meant to be.
# The depths must agree exactly. The values will not, because its subsample is unseeded, so
# the check is that two independent draws at the same depth agree on richness to within
# subsampling noise, not that they are equal.
x_rarefied_qiime <- ap_import(ap_read_qza(inp[["rarefied"]]), meta, tree = tree,
                              taxonomy = taxonomy)
q0_r <- ap_alpha(x_rarefied, metrics = "q0", rarefy = FALSE)$values
q0_q <- ap_alpha(x_rarefied_qiime, metrics = "q0", rarefy = FALSE)$values
shared <- intersect(q0_r$sample, q0_q$sample)
cross <- data.frame(
  quantity = c("samples_r", "samples_qiime", "samples_shared", "depth_r", "depth_qiime",
               "q0_mean_abs_diff", "q0_max_abs_diff", "q0_spearman"),
  value = c(
    ncol(x_rarefied), ncol(x_rarefied_qiime), length(shared),
    unique(colSums(SummarizedExperiment::assay(x_rarefied, "counts"))),
    unique(colSums(SummarizedExperiment::assay(x_rarefied_qiime, "counts"))),
    mean(abs(q0_r$value[match(shared, q0_r$sample)] - q0_q$value[match(shared, q0_q$sample)])),
    max(abs(q0_r$value[match(shared, q0_r$sample)] - q0_q$value[match(shared, q0_q$sample)])),
    stats::cor(q0_r$value[match(shared, q0_r$sample)],
               q0_q$value[match(shared, q0_q$sample)], method = "spearman")
  )
)
print(cross)
write_tsv(cross, "rarefaction_crosscheck")

# --- design -----------------------------------------------------------------------------------

section("metadata scan")
scan <- ap_scan_metadata(x, group = an$group)
print(scan)
write_tsv(scan$variables, "metadata_variables")
if (!is.null(scan$confounders)) write_tsv(scan$confounders, "metadata_confounders")
write_tsv(ap_depth_candidates(x), "depth_candidates")

# --- alpha diversity ----------------------------------------------------------------------------

section("alpha diversity")
alpha <- ap_alpha(x_rarefied, metrics = c("q0", "q1", "q2", "evenness", "faith_pd"),
                  rarefy = FALSE)
alpha_test <- ap_alpha_test(alpha, an$group, seed = an$seed)
print(alpha_test)
write_tsv(alpha$values, "alpha_values")
write_tsv(alpha_test$results, "alpha_tests")
try_plot(ap_plot_alpha(alpha, an$group, test = alpha_test), "alpha", width = 9, height = 6)

# Repeated rarefaction of the unrarefied table, for comparison with the single subsample
# QIIME 2 draws. Reported alongside, not substituted.
alpha_repeated <- ap_alpha(x, metrics = c("q0", "q1", "q2", "evenness", "faith_pd"),
                           rarefy = TRUE, depth = depth,
                           n_iter = an$alpha_rarefy_iterations, seed = an$seed)
write_tsv(alpha_repeated$values, "alpha_values_repeated_rarefaction")

# --- beta diversity -----------------------------------------------------------------------------

section("beta diversity")
beta <- ap_beta(x_rarefied, metrics = c("bray_curtis", "jaccard", "unweighted_unifrac",
                                        "weighted_unifrac"), rarefy = FALSE)
permanova <- ap_permanova(beta, an$group, permutations = an$permutations, seed = an$seed)
print(permanova)
write_tsv(permanova$results, "permanova")
write_tsv(permanova$dispersion, "permanova_dispersion")
write_tsv(permanova$interpretation, "permanova_interpretation")
for (m in beta$metrics) {
  ord <- ap_ordinate(beta, metric = m, seed = an$seed)
  try_plot(ap_plot_ordination(ord, group = an$group, permanova = permanova),
           paste0("ordination_", m))
  try_plot(ap_plot_dispersion(permanova, metric = m, term = an$group), paste0("dispersion_", m))
}

# --- composition --------------------------------------------------------------------------------

section("composition")
write_tsv(ap_top_taxa(x, n = 15L, group = an$group, rank = "genus"), "top_genera")
try_plot(ap_plot_taxa_bar(x, rank = "phylum", n = 10L, group = an$group, mode = "group"),
         "taxa_phylum_bars")
try_plot(ap_plot_taxa_bar(x, rank = "genus", n = 15L, group = an$group, mode = "group"),
         "taxa_genus_bars", width = 9)
try_plot(ap_plot_taxa_heatmap(x, rank = "genus", n = 25L, group = an$group),
         "taxa_genus_heatmap", width = 10, height = 6)

# --- differential abundance ---------------------------------------------------------------------

section("differential abundance")
x_da <- x
if (length(an$da_levels) > 0L) {
  g <- SummarizedExperiment::colData(x)[[an$group]]
  x_da <- x[, !is.na(g) & g %in% an$da_levels]
  cat("differential abundance on", ncol(x_da), "samples with", an$group, "in",
      paste(an$da_levels, collapse = ", "), "\n")
}
da_args <- list(x_da, group = an$group, reference = an$da_reference, seed = an$seed)
if (length(an$da_methods) > 0L) da_args$method <- unlist(an$da_methods)
da <- do.call(ap_da, da_args)
print(da)
write_tsv(da$results, "da_results")
concordance <- ap_da_concordance(da)
print(concordance)
write_tsv(concordance$features, "da_concordance")
try_plot(ap_plot_volcano(da), "da_volcano", width = 9, height = 6)
try_plot(ap_plot_concordance(concordance), "da_concordance", width = 8, height = 8)
try_plot(ap_plot_da_effects(concordance), "da_effects", width = 8, height = 6)

# --- screen and confirmatory model ---------------------------------------------------------------

section("screen")
screen <- ap_screen(alpha = alpha, beta = beta, permutations = an$permutations,
                    n_resample = an$n_resample, seed = an$seed)
print(screen, n = 30L)
write_tsv(screen$results, "screen_results")
write_tsv(screen$by_variable, "screen_by_variable")
if (!is.null(screen$skipped)) write_tsv(screen$skipped, "screen_skipped")
try_plot(ap_plot_screen(screen, n = 25L), "screen", width = 6.5, height = 5.5)

explains <- NULL
if (length(an$explains_terms) >= 2L) {
  section("explains")
  explains <- ap_explains(beta = beta, alpha = alpha, terms = unlist(an$explains_terms),
                          permutations = an$permutations, seed = an$seed)
  print(explains)
  write_tsv(explains$beta$terms, "explains_beta_terms")
  write_tsv(explains$beta$model, "explains_beta_model")
  write_tsv(explains$alpha$terms, "explains_alpha_terms")
  write_tsv(explains$alpha$model, "explains_alpha_model")
}

# --- normalization sensitivity -------------------------------------------------------------------

section("normalization sensitivity")
normalization <- ap_normalization_sensitivity(
  x, terms = an$group, methods = unlist(an$normalizations), depth = depth,
  permutations = an$permutations, seed = an$seed
)
print(normalization)
write_tsv(normalization$results, "normalization_results")
write_tsv(normalization$agreement, "normalization_agreement")

# --- results and provenance --------------------------------------------------------------------

section("provenance")
saveRDS(list(scan = scan, depth = depth, alpha = alpha, alpha_test = alpha_test,
             alpha_repeated = alpha_repeated, beta = beta, permanova = permanova,
             da = da, concordance = concordance, screen = screen, explains = explains,
             normalization = normalization, config = cfg),
        snakemake@output[["results"]])

sha256 <- function(path) {
  p <- normalizePath(path)
  sub(" .*$", "", system2("sha256sum", shQuote(p), stdout = TRUE))
}
# snakemake@input holds every input twice, once positionally and once by name, so the
# unnamed half would write six sha256 rows labelled only "sha256:" with no file behind them.
input_files <- unlist(inp)
input_files <- input_files[nzchar(names(input_files))]
run_info <- data.frame(
  key = c("started", "finished", "r_version", "amplipub_version", "seed", "permutations",
          "n_resample", "rarefaction_depth", paste0("sha256:", names(input_files))),
  value = c(format(started, "%Y-%m-%d %H:%M:%S %Z"), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            R.version.string, as.character(utils::packageVersion("AmpliPub")),
            an$seed, an$permutations, an$n_resample, depth,
            vapply(input_files, sha256, character(1))),
  stringsAsFactors = FALSE
)
utils::write.table(run_info, snakemake@output[["run_info"]], sep = "\t", quote = FALSE,
                   row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(dir_prov, "sessionInfo.txt"))
yaml::write_yaml(cfg, file.path(dir_prov, "config_resolved.yaml"))
cat("done in", format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")
