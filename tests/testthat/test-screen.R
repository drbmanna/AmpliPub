# The screen exists to stop two things: ranking by p, and multiple testing that
# nobody counts. Each test below is built so the right answer is known first.

test_that("effect size sets the order, not the p-value", {
  # `a` is a modest shift measured on all 200 samples, so its p is tiny. `b` is a
  # perfect separation measured on only 16 samples, so its effect is far larger
  # and its p is not. Ranking by p would put `a` first.
  set.seed(101)
  n <- 200L
  meta <- data.frame(a = rep(c("lo", "hi"), each = n / 2),
                     row.names = sprintf("S%03d", seq_len(n)), stringsAsFactors = FALSE)
  meta$b <- NA_character_
  meta$b[1:8] <- "x"
  meta$b[101:108] <- "y"
  y <- stats::rnorm(n) + ifelse(meta$a == "hi", 1, 0)
  y[101:108] <- y[101:108] + 10

  s <- ap_screen(alpha = ap_fixture_alpha(y, meta), top_k = 1L, n_resample = 10L)
  ra <- s$results[s$results$variable == "a", ]
  rb <- s$results[s$results$variable == "b", ]
  expect_lt(ra$p, rb$p)
  expect_gt(rb$effect, ra$effect)
  expect_equal(s$results$variable[1], "b")
})

test_that("one BH correction runs across the whole screen and the test count is recorded", {
  f <- ap_fixture_screen_inputs()
  s <- ap_screen(alpha = f$alpha, beta = f$beta, permutations = 99L, n_resample = 20L)
  # group, batch_run, sex, age across q0, q1 and bray_curtis.
  expect_equal(s$n_tests, 12L)
  expect_equal(nrow(s$results), s$n_tests)
  expect_equal(s$results$q, stats::p.adjust(s$results$p, method = "BH"))
  expect_equal(s$expected_false_positives, 0.05 * 12)
  expect_equal(s$results$effect_adj, sort(s$results$effect_adj, decreasing = TRUE))
})

test_that("printing states the test count and that the screen is hypothesis-generating", {
  f <- ap_fixture_screen_inputs()
  s <- ap_screen(alpha = f$alpha, beta = f$beta, permutations = 99L, n_resample = 20L)
  stdout_txt <- paste(utils::capture.output(
    msg_txt <- paste(utils::capture.output(print(s), type = "message"), collapse = " ")
  ), collapse = " ")
  expect_match(msg_txt, "hypothesis-generating", fixed = TRUE)
  expect_match(msg_txt, "12 tests", fixed = TRUE)
  expect_match(stdout_txt, "stable")
})

test_that("a strongly planted variable is stable, and stability shares add up to top_k", {
  set.seed(202)
  n <- 120L
  meta <- data.frame(strong = rep(c("a", "b"), each = n / 2),
                     row.names = sprintf("S%03d", seq_len(n)), stringsAsFactors = FALSE)
  for (i in 1:6) meta[[paste0("noise", i)]] <- sample(rep(c("u", "v"), each = n / 2))
  y <- stats::rnorm(n) + ifelse(meta$strong == "b", 3, 0)

  s <- ap_screen(alpha = ap_fixture_alpha(y, meta), top_k = 3L, n_resample = 50L, seed = 5L)
  expect_equal(s$results$variable[1], "strong")
  expect_gte(s$results$stability[s$results$variable == "strong"], 0.9)
  # Exactly top_k rows are in the top k of every subsample when effects do not
  # tie, so the shares sum to top_k. A bookkeeping error breaks this.
  expect_lt(abs(sum(s$results$stability) - 3), 0.1)
})

test_that("subsamples are drawn without replacement at the stated fraction", {
  f <- ap_fixture_screen_inputs()
  s <- ap_screen(alpha = f$alpha, n_resample = 20L)
  expect_length(s$resamples, 20L)
  expect_true(all(vapply(s$resamples, function(r) length(r) == 19L, logical(1))))
  expect_false(any(vapply(s$resamples, anyDuplicated, integer(1)) > 0L))
})

test_that("settings that would make stability meaningless are refused", {
  f <- ap_fixture_screen_inputs()
  expect_error(ap_screen(alpha = f$alpha, fraction = 1), "strictly between")
  expect_error(ap_screen(alpha = f$alpha, n_resample = 5L), "at least 10")
  expect_error(ap_screen(alpha = f$alpha, top_k = 50L, n_resample = 10L),
               "smaller than the number of tests")
  expect_error(ap_screen(), "at least one")
})

test_that("the same seed gives the same screen", {
  f <- ap_fixture_screen_inputs()
  s1 <- ap_screen(alpha = f$alpha, beta = f$beta, permutations = 99L, n_resample = 20L, seed = 9L)
  s2 <- ap_screen(alpha = f$alpha, beta = f$beta, permutations = 99L, n_resample = 20L, seed = 9L)
  expect_identical(s1$results, s2$results)
  expect_identical(s1$resamples, s2$resamples)
})

test_that("identifiers, constants, repeated measures and tiny levels are skipped with a reason", {
  set.seed(303)
  n <- 100L
  meta <- data.frame(
    grp = rep(c("a", "b"), each = n / 2),
    age = round(stats::runif(n, 20, 80)),
    sample_code = sprintf("id%03d", seq_len(n)),
    fixed = "same",
    subject = rep(sprintf("p%02d", 1:50), each = 2),
    rare = c(rep("common", n - 2), "odd", "odd"),
    row.names = sprintf("S%03d", seq_len(n)), stringsAsFactors = FALSE
  )
  obj <- ap_fixture_alpha(stats::rnorm(n), meta)
  s <- ap_screen(alpha = obj, top_k = 1L, n_resample = 10L)

  sk <- s$skipped
  expect_setequal(sk$variable, c("sample_code", "fixed", "subject", "rare"))
  expect_match(sk$reason[sk$variable == "sample_code"], "identifier")
  expect_match(sk$reason[sk$variable == "fixed"], "constant")
  expect_match(sk$reason[sk$variable == "subject"], "mixed model")
  expect_match(sk$reason[sk$variable == "rare"], "fewer than 3")
  expect_setequal(s$results$variable, c("grp", "age"))

  expect_error(ap_screen(alpha = obj, variables = "nope"), "not in the sample metadata")
  expect_error(ap_screen(alpha = obj, variables = "fixed"), "No variable is testable")
})

test_that("beta rows carry the dispersion verdict and a confound shows up as an equal effect", {
  f <- ap_fixture_screen_inputs()
  s <- ap_screen(beta = f$beta, permutations = 99L, n_resample = 20L, top_k = 2L)
  b <- s$results
  expect_true(all(!is.na(b$verdict)))
  expect_false(is.na(b$dispersion_p[b$variable == "group"]))
  # batch_run is identical to group in the fixture, so nothing can tell them apart.
  expect_equal(b$effect[b$variable == "group"], b$effect[b$variable == "batch_run"])
})

test_that("the fast R2 used for stability equals adonis2 exactly", {
  f <- ap_fixture_screen_inputs()
  d <- f$beta$distances$bray_curtis
  m <- as.matrix(d)
  meta <- f$beta$metadata[rownames(m), ]

  g <- factor(meta$group)
  ref <- vegan::adonis2(d ~ g, data = data.frame(g = g), permutations = 0)[1, "R2"]
  expect_equal(AmpliPub:::ap_screen_r2(g, m = m), ref, tolerance = 1e-10)

  age <- meta$age
  ref <- vegan::adonis2(d ~ age, data = data.frame(age = age), permutations = 0)[1, "R2"]
  expect_equal(AmpliPub:::ap_screen_r2(age, m = m), ref, tolerance = 1e-10)

  # A Gower matrix built on all samples is wrong once a sample is missing. The
  # function must notice and rebuild rather than reuse it.
  g[1:3] <- NA
  keep <- !is.na(g)
  dd <- stats::as.dist(m[keep, keep])
  ref <- vegan::adonis2(dd ~ g, data = data.frame(g = g[keep]), permutations = 0)[1, "R2"]
  stale <- AmpliPub:::ap_gower(m)
  expect_equal(AmpliPub:::ap_screen_r2(g, G = stale, m = m), ref, tolerance = 1e-10)
})

test_that("the screen plot builds", {
  f <- ap_fixture_screen_inputs()
  s <- ap_screen(alpha = f$alpha, beta = f$beta, permutations = 99L, n_resample = 20L)
  p <- ap_plot_screen(s, n = 8L)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_error(ap_plot_screen(list()), "must come from")
})

test_that("a null variable with many levels does not outrank a real two-level effect", {
  # Found on Baxter: raw R2 has a null expectation near df/(n - 1). Here a
  # 250-level variable with no effect sits near 0.25 raw, above a planted
  # two-group effect near 0.14. Ranking on raw R2 puts noise first; the
  # adjusted scale must not.
  #
  # Sized so a single draw decides it. At n = 200 with 50 levels the null's
  # adjusted R2 has sd 0.055 over label permutations (checked), so one draw
  # landed 2 sd high and outranked the real effect. At n = 1000 both
  # comparisons below sit about 5 sd from failing.
  set.seed(404)
  n <- 1000L
  real <- rep(c("a", "b"), each = n / 2)
  y <- stats::rnorm(n) + ifelse(real == "b", 0.8, 0)
  meta <- data.frame(real = real,
                     many = sample(rep(sprintf("L%03d", 1:250), each = 4)),
                     row.names = sprintf("S%04d", seq_len(n)), stringsAsFactors = FALSE)
  d <- stats::dist(matrix(y, ncol = 1, dimnames = list(rownames(meta), NULL)))
  beta <- structure(list(distances = list(euclidean = d), metrics = "euclidean",
                         rarefied = FALSE, depth = NA_real_, seed = NA_integer_,
                         dropped = character(0), pseudocount = NA_real_,
                         n_samples = n, metadata = meta),
                    class = "ap_beta")

  # No assertion here reads a p-value, so permutations are kept minimal; at
  # n = 1000 with 250 levels each one is expensive.
  s <- ap_screen(beta = beta, permutations = 19L, n_resample = 10L, top_k = 1L,
                 max_levels = 300L)
  r <- s$results
  # The fixture does what it claims: noise wins on the raw scale...
  expect_gt(r$effect[r$variable == "many"], r$effect[r$variable == "real"])
  # ...and the ranking does not follow it.
  expect_equal(r$variable[1], "real")
  expect_equal(r$df[r$variable == "many"], 249)
  expect_equal(r$effect_adj, 1 - (1 - r$effect) * (r$n - 1) / (r$n - r$df - 1))
})
