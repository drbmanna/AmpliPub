# ap_explains() partitions variance between terms. The fixtures plant overlap
# and independence so the right partition is known before fitting.

# `near_group` agrees with `group` on 20 of 24 samples: strongly overlapping but
# not collinear, so both marginal terms exist and most variance is shared.
ap_fixture_explains <- function() {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  cd <- SummarizedExperiment::colData(x)
  near <- cd$group
  flip <- c(1, 2, 13, 14)
  near[flip] <- ifelse(near[flip] == "a", "b", "a")
  cd$near_group <- near
  SummarizedExperiment::colData(x) <- cd
  x
}

test_that("a noisy proxy keeps almost nothing of its own once the real term is in", {
  # The signal is planted on `group`; `near_group` only carries it by agreeing
  # with group on 20 of 24 samples. So near_group's marginal R2 should be close
  # to zero and what it appears to explain on its own should be shared.
  x <- ap_fixture_explains()
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  e <- ap_explains(beta = b, terms = c("near_group", "group"), permutations = 99L)
  t <- e$beta$terms
  expect_lt(t$R2[t$term == "near_group"], 0.01)
  expect_gt(e$beta$model$shared_R2, 10 * t$R2[t$term == "near_group"])

  # Two-term identity: a term fitted alone explains its marginal R2 plus the
  # shared part. Checked against a separate single-term PERMANOVA.
  alone <- ap_permanova(b, "near_group", permutations = 99L)$results$R2
  expect_equal(alone, t$R2[t$term == "near_group"] + e$beta$model$shared_R2,
               tolerance = 1e-10)
})

test_that("the model R2 agrees with adonis2 fitting all terms at once", {
  x <- ap_fixture_explains()
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  e <- ap_explains(beta = b, terms = c("near_group", "group"), permutations = 99L)

  d <- b$distances$bray_curtis
  md <- b$metadata[attr(d, "Labels"), c("near_group", "group")]
  md$near_group <- factor(md$near_group)
  md$group <- factor(md$group)
  ref <- vegan::adonis2(d ~ near_group + group, data = md, permutations = 0)[1, "R2"]
  expect_equal(e$beta$model$model_R2, ref, tolerance = 1e-10)
})

test_that("term order does not change any term's R2 or the shared part", {
  x <- ap_fixture_explains()
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  e1 <- ap_explains(beta = b, terms = c("near_group", "group"), permutations = 99L)
  e2 <- ap_explains(beta = b, terms = c("group", "near_group"), permutations = 99L)
  expect_equal(e2$beta$terms$R2[match(e1$beta$terms$term, e2$beta$terms$term)],
               e1$beta$terms$R2)
  expect_equal(e2$beta$model$shared_R2, e1$beta$model$shared_R2)
})

test_that("independent terms in a balanced design share nothing", {
  # In the fixture `sex` alternates within each group, six of each, so group and
  # sex are orthogonal and their sums of squares add exactly.
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  a <- suppressMessages(ap_alpha(x, metrics = "q1", rarefy = FALSE))
  e <- ap_explains(beta = b, alpha = a, terms = c("group", "sex"), permutations = 99L)
  expect_lt(abs(e$beta$model$shared_R2), 1e-8)
  expect_lt(abs(e$alpha$model$shared_R2), 1e-8)
})

test_that("alpha marginal R2 equals the drop in residual sum of squares, computed separately", {
  x <- ap_fixture_explains()
  a <- suppressMessages(ap_alpha(x, metrics = "q1", rarefy = FALSE))
  e <- ap_explains(alpha = a, terms = c("near_group", "group"))

  df <- data.frame(value = a$values$value,
                   a$metadata[a$values$sample_id, c("near_group", "group")])
  full <- stats::lm(value ~ near_group + group, data = df)
  reduced <- stats::lm(value ~ near_group, data = df)
  sst <- sum((df$value - mean(df$value))^2)
  t <- e$alpha$terms
  expect_equal(t$R2[t$term == "group"],
               (stats::deviance(reduced) - stats::deviance(full)) / sst, tolerance = 1e-10)
  expect_equal(e$alpha$model$model_R2, summary(full)$r.squared)
})

test_that("perfectly collinear terms are refused rather than reported as zero", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  a <- suppressMessages(ap_alpha(x, metrics = "q1", rarefy = FALSE))
  expect_error(ap_explains(beta = b, terms = c("batch_run", "group"), permutations = 99L),
               "perfectly collinear")
  expect_error(ap_explains(alpha = a, terms = c("batch_run", "group")), "perfectly collinear")
})

test_that("bad input is refused", {
  x <- ap_fixture_object(tree = FALSE, taxonomy = FALSE)
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  expect_error(ap_explains(terms = c("group", "sex")), "at least one")
  expect_error(ap_explains(beta = b, terms = "group"), "at least two")
  expect_error(ap_explains(beta = b, terms = c("group", "nope")), "not in the sample metadata")
})

test_that("printing reports the shared variance", {
  x <- ap_fixture_explains()
  b <- ap_beta(x, metrics = "bray_curtis", rarefy = FALSE)
  e <- ap_explains(beta = b, terms = ~ near_group + group, permutations = 99L)
  msg_txt <- paste(utils::capture.output(
    invisible(utils::capture.output(print(e))), type = "message"), collapse = " ")
  expect_match(msg_txt, "shared", ignore.case = TRUE)
  expect_match(msg_txt, "Confirmatory model", fixed = TRUE)
})
