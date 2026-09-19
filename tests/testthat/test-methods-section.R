# The methods section states what was done to this dataset, so unlike the
# reader's guide it is expected to vary with the data. What must never happen is
# a number, version or citation appearing that did not come from the run.

test_that("no function is defined twice across R/", {
  # Written after ap_metric_label was defined in both theme.R and
  # methods_section.R. The duplicate did not error: R simply kept whichever file
  # loaded last, and the other definition vanished. Here the plot labels won and
  # the methods prose silently printed raw metric names. Had the file been named
  # differently it would have gone the other way and broken every plot label.
  r_dir <- testthat::test_path("..", "..", "R")
  if (!dir.exists(r_dir)) skip("source R/ not present (installed package)")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  if (!length(files)) skip("no R sources found")
  defs <- unlist(lapply(files, function(f) {
    ln <- readLines(f, warn = FALSE)
    sub(" <- function.*$", "", grep("^[a-zA-Z_.][a-zA-Z0-9_.]* <- function", ln, value = TRUE))
  }))
  dup <- unique(defs[duplicated(defs)])
  expect_equal(dup, character(0))
})

# --- version parsing ----------------------------------------------------------

test_that("conda versions are read from the environment export", {
  d <- withr::local_tempdir()
  writeLines(c("name: test", "dependencies:", "  - cutadapt=5.1",
               "  - mafft=7.526", "  - q2-dada2=2025.7.0", "  - fasttree=2.1.11=h7b50bb2_0"),
             file.path(d, "env_test.yml"))
  v <- ap_env_versions(d)
  expect_equal(v$cutadapt, "5.1")
  expect_equal(v$mafft, "7.526")
  expect_equal(v$fasttree, "2.1.11")
  expect_equal(v[["qiime2-version"]], "2025.7")
})

test_that("a missing or unreadable provenance directory yields no versions, not guesses", {
  expect_equal(ap_env_versions(NULL), list())
  expect_equal(ap_env_versions(file.path(tempdir(), "does-not-exist")), list())
})

# --- key/value files ----------------------------------------------------------

test_that("headerless key/value summaries are read without losing the first row", {
  d <- withr::local_tempdir()
  dir.create(file.path(d, "q2", "quality"), recursive = TRUE)
  writeLines(c("trunc_len_f\t251", "trunc_len_r\t219", "expected_overlap\t217"),
             file.path(d, "q2", "quality", "trunc_len.tsv"))
  tl <- ap_read_run_tsv(d, "q2/quality/trunc_len.tsv", header = FALSE)
  expect_equal(ap_kv(tl, "trunc_len_f"), "251")
  expect_equal(ap_kv(tl, "expected_overlap"), "217")
})

test_that("a missing run file returns NULL rather than failing the whole section", {
  expect_null(ap_read_run_tsv(NULL, "q2/quality/trunc_len.tsv"))
  expect_null(ap_read_run_tsv(tempdir(), "q2/nope/missing.tsv"))
  expect_true(is.na(ap_kv(NULL, "anything")))
})

# --- citations ----------------------------------------------------------------

test_that("inline citations come from the package's own CITATION", {
  key <- ap_cite_inline("vegan")
  expect_true(grepl("[0-9]{4}$", key))
  expect_false(grepl("^vegan", key))
})

test_that("a package with no usable citation falls back to name and version, never invention", {
  key <- ap_cite_inline("this.package.does.not.exist")
  expect_equal(key, "this.package.does.not.exist")
})

test_that("reference entries are stripped of formatting artefacts but not of content", {
  # The asterisks are followed by a line break here, not a space, which is how
  # ALDEx2's own CITATION is laid out.
  s <- ap_clean_reference(c('Smith A (2020). "***', 'A Title." _*** J Name_.',
                            "<***%20http://doi.org/10.1/x>"))
  expect_false(grepl('" A Title', s, fixed = TRUE))
  expect_false(grepl("*", s, fixed = TRUE))
  expect_true(grepl('"A Title."', s, fixed = TRUE))
  expect_true(grepl("http://doi.org/10.1/x", s, fixed = TRUE))
  expect_true(grepl("Smith A (2020)", s, fixed = TRUE))
  # The gap after the closing quote must survive, or the title runs into the
  # journal name.
  expect_true(grepl('Title." _J Name_', s, fixed = TRUE))
})

# --- prose labels -------------------------------------------------------------

test_that("distances and ordinations are named as a reader expects, not as object keys", {
  expect_equal(ap_distance_label("bray_curtis"), "Bray-Curtis dissimilarity")
  expect_equal(ap_distance_label("unweighted_unifrac"), "unweighted UniFrac")
  expect_equal(ap_ord_label("pcoa"), "principal coordinates analysis (PCoA)")
  expect_equal(ap_da_label("ancombc2"), "ANCOM-BC2")
  expect_equal(ap_norm_label("clr"), "a centred log-ratio transform with Aitchison distance")
})

test_that("an unrecognised key degrades to itself rather than erroring", {
  expect_equal(ap_distance_label("some_new_metric"), "some new metric")
  expect_equal(ap_da_label("newmethod"), "newmethod")
})

test_that("ap_and joins one, two and many items correctly", {
  expect_equal(ap_and(character(0)), "")
  expect_equal(ap_and("a"), "a")
  expect_equal(ap_and(c("a", "b")), "a and b")
  expect_equal(ap_and(c("a", "b", "c")), "a, b and c")
})

# --- assembly -----------------------------------------------------------------

ap_methods_res <- function(...) {
  utils::modifyList(
    list(config = list(analysis = list(seed = 1L, permutations = 999L)),
         depth = NULL, alpha_test = NULL, beta = NULL, permanova = NULL,
         ordination_diagnostics = NULL, da = NULL, concordance = NULL,
         screen = NULL, explains = NULL, normalization = NULL),
    list(...))
}

test_that("only the analyses that ran are described", {
  m <- paste(ap_methods_section(ap_methods_res()), collapse = "\n")
  expect_true(grepl("# Methods", m, fixed = TRUE))
  expect_false(grepl("PERMANOVA", m, fixed = TRUE))
  expect_false(grepl("Differential abundance was assessed", m, fixed = TRUE))
  expect_false(grepl("exploratory screen", m, fixed = TRUE))
})

test_that("numbers that cannot be read are omitted rather than invented", {
  # No run_dir, so no truncation lengths, read counts or filtered table sizes.
  res <- ap_methods_res(config = list(
    analysis = list(seed = 1L),
    quality = list(n = 10000L, min_q = 30L, min_overlap = 12L,
                   amplicon_len = 253L, margin = 20L)))
  m <- paste(ap_methods_section(res), collapse = "\n")
  expect_true(grepl("subsampled reads", m, fixed = TRUE))
  expect_false(grepl("were truncated at", m, fixed = TRUE))
  expect_false(grepl("input read pairs", m, fixed = TRUE))
})

test_that("a cited package appears in the reference list exactly once", {
  res <- ap_methods_res(alpha_test = list(x = 1), permanova = list(x = 1),
                        beta = list(metrics = "bray_curtis"))
  m <- ap_methods_section(res)
  refs <- grep("^- ", m, value = TRUE)
  expect_equal(sum(grepl("vegan", refs)), 1L)
})

test_that("the draft is labelled as a draft, because it is not submission-ready", {
  m <- paste(ap_methods_section(ap_methods_res()), collapse = "\n")
  expect_true(grepl("Draft generated", m, fixed = TRUE))
  expect_true(grepl("Check every number", m, fixed = TRUE))
})

test_that("the preprocessing summary names each analysis's table and normalization", {
  res <- list(
    depth = 5000L,
    alpha = list(), alpha_repeated = list(rarefied = TRUE, n_iter = 10L, depth = 5000L, seed = 3L),
    beta = list(metrics = c("bray_curtis", "aitchison"), pseudocount = 0.5),
    da = list(prv_cut = 0.1, n_features = 50L, n_samples = 20L,
              methods = c("ancombc2", "maaslin2")),
    normalization = list(methods = c("tss", "rarefy"), permutations = 99L),
    config = list(analysis = list(seed = 3L, da_levels = list("a", "b")),
                  filter = list(include = "Bacteria", exclude = "mitochondria,chloroplast",
                                min_samples_fraction = 0.05))
  )
  s <- ap_preprocessing_summary(res)
  expect_named(s, c("analysis", "input", "filter", "normalization"))
  expect_match(s$normalization[s$analysis == "Alpha diversity and its tests"], "5,000 reads in R, seed 3")
  beta <- s[grepl("^Beta diversity", s$analysis), ]
  expect_match(beta$normalization, "Aitchison: zeros replaced by 0.5")
  da <- s[grepl("^Differential abundance", s$analysis), ]
  expect_match(da$input, "raw counts")
  expect_match(da$filter, "10% of those samples")
  expect_match(da$normalization, "pseudocount sensitivity")
  expect_no_match(da$normalization, "ALDEx2")
  expect_match(s$filter[1], "mitochondria and chloroplast sequences removed")
  expect_match(s$normalization[s$analysis == "Normalization sensitivity"], "TSS and rarefying")
})
