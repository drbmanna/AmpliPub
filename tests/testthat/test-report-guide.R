# The guide must describe every analysis that ran, describe none that did not,
# and say the same thing regardless of what the numbers turned out to be.

ap_guide_res <- function(..., scan_extras = list()) {
  scan <- utils::modifyList(
    list(variables = data.frame(variable = "dx", role = "categorical",
                                stringsAsFactors = FALSE),
         batch_candidates = character(0),
         repeated_measures = character(0),
         confounders = NULL,
         group = "dx", n_samples = 100L),
    scan_extras)
  utils::modifyList(list(scan = scan, depth = NULL, alpha_test = NULL,
                         permanova = NULL, ordination_diagnostics = NULL,
                         da = NULL, concordance = NULL, screen = NULL,
                         explains = NULL, normalization = NULL),
                    list(...))
}

guide_text <- function(res) paste(ap_report_guide(res), collapse = "\n")

test_that("the frame, study design, composition and limits are always present", {
  g <- guide_text(ap_guide_res())
  expect_true(grepl("How to read this report", g, fixed = TRUE))
  expect_true(grepl("## Study design", g, fixed = TRUE))
  expect_true(grepl("## Composition", g, fixed = TRUE))
  expect_true(grepl("## What this report does not do", g, fixed = TRUE))
})

test_that("an analysis that did not run is not described", {
  g <- guide_text(ap_guide_res())
  expect_false(grepl("## Alpha diversity", g, fixed = TRUE))
  expect_false(grepl("## Beta diversity", g, fixed = TRUE))
  expect_false(grepl("## Ordination plots", g, fixed = TRUE))
  expect_false(grepl("## Differential abundance", g, fixed = TRUE))
  expect_false(grepl("## Exploratory screen", g, fixed = TRUE))
  expect_false(grepl("## Confirmatory model", g, fixed = TRUE))
  expect_false(grepl("## Normalization sensitivity", g, fixed = TRUE))
  expect_false(grepl("## Sequencing depth", g, fixed = TRUE))
})

test_that("each analysis that ran adds exactly its own section", {
  expect_true(grepl("## Alpha diversity",
                    guide_text(ap_guide_res(alpha_test = list(x = 1))), fixed = TRUE))
  expect_true(grepl("## Beta diversity",
                    guide_text(ap_guide_res(permanova = list(x = 1))), fixed = TRUE))
  expect_true(grepl("## Ordination plots",
                    guide_text(ap_guide_res(ordination_diagnostics = data.frame(a = 1))),
                    fixed = TRUE))
  expect_true(grepl("## Differential abundance",
                    guide_text(ap_guide_res(da = list(x = 1))), fixed = TRUE))
  expect_true(grepl("## Exploratory screen",
                    guide_text(ap_guide_res(screen = list(x = 1))), fixed = TRUE))
  expect_true(grepl("## Confirmatory model",
                    guide_text(ap_guide_res(explains = list(x = 1))), fixed = TRUE))
  expect_true(grepl("## Normalization sensitivity",
                    guide_text(ap_guide_res(normalization = list(x = 1))), fixed = TRUE))
  expect_true(grepl("## Sequencing depth",
                    guide_text(ap_guide_res(depth = 10000)), fixed = TRUE))
})

test_that("the dispersion paragraph appears only when dispersion was tested", {
  with_disp <- guide_text(ap_guide_res(permanova = list(dispersion = data.frame(a = 1))))
  without <- guide_text(ap_guide_res(permanova = list(results = data.frame(a = 1))))
  expect_true(grepl("betadisper", with_disp, fixed = TRUE))
  expect_false(grepl("betadisper", without, fixed = TRUE))
})

test_that("design warnings appear only when the scan raised them", {
  # "Candidate confounders" also appears in the always-present limits block, so
  # each paragraph is matched on a phrase unique to it.
  plain <- guide_text(ap_guide_res())
  expect_false(grepl("Batch and technical proxies", plain, fixed = TRUE))
  expect_false(grepl("are not independent", plain, fixed = TRUE))
  expect_false(grepl("chi-squared test of independence", plain, fixed = TRUE))

  batch <- guide_text(ap_guide_res(scan_extras = list(batch_candidates = "run")))
  expect_true(grepl("Batch and technical proxies", batch, fixed = TRUE))

  rep <- guide_text(ap_guide_res(scan_extras = list(repeated_measures = "subject")))
  expect_true(grepl("are not independent", rep, fixed = TRUE))

  conf <- guide_text(ap_guide_res(scan_extras = list(confounders = data.frame(a = 1))))
  expect_true(grepl("chi-squared test of independence", conf, fixed = TRUE))
})

test_that("the Hill number naming is spelled out, because users look for Shannon and Simpson", {
  g <- guide_text(ap_guide_res(alpha_test = list(x = 1)))
  expect_true(grepl("This is Shannon diversity", g, fixed = TRUE))
  expect_true(grepl("This is Simpson diversity", g, fixed = TRUE))
  expect_true(grepl("`q1`", g, fixed = TRUE))
  expect_true(grepl("`q2`", g, fixed = TRUE))
})

test_that("the screen is described as hypothesis-generating, not as findings", {
  g <- guide_text(ap_guide_res(screen = list(x = 1)))
  expect_true(grepl("candidates, not findings", g, fixed = TRUE))
  expect_true(grepl("garden of forking paths", g, fixed = TRUE))
})

test_that("the guide does not change with the data, only with which analyses ran", {
  # Same structure, wildly different contents. The guide is about methods, so it
  # must be byte-identical. If this fails, something data-specific leaked in.
  a <- ap_guide_res(depth = 10000, alpha_test = list(results = data.frame(p = 0.001)),
                    permanova = list(dispersion = data.frame(p = 0.02)),
                    ordination_diagnostics = data.frame(stress = 0.05),
                    da = list(n = 500), concordance = list(k = 4),
                    screen = list(hits = 90), explains = list(r2 = 0.9),
                    normalization = list(m = 4))
  b <- ap_guide_res(depth = 2500, alpha_test = list(results = data.frame(p = 0.9)),
                    permanova = list(dispersion = data.frame(p = 0.77)),
                    ordination_diagnostics = data.frame(stress = 0.31),
                    da = list(n = 0), concordance = list(k = 2),
                    screen = list(hits = 0), explains = list(r2 = 0.01),
                    normalization = list(m = 4))
  expect_identical(ap_report_guide(a), ap_report_guide(b))
})

test_that("no sample count, p-value or metric value leaks into the guide", {
  g <- guide_text(ap_guide_res(depth = 12345, alpha_test = list(x = 1),
                               permanova = list(x = 1), screen = list(x = 1)))
  expect_false(grepl("12345", g, fixed = TRUE))
  expect_false(grepl("dx", g, fixed = TRUE))
})
