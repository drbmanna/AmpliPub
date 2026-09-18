# A reader's guide, emitted at the end of the report.
#
# The per-result interpretation text elsewhere in this package says what a
# particular number means. This says what the analysis was, what its columns
# are, and what a reader may not claim from it. It is the part a user who did
# not design the study needs in order to read their own output.
#
# Every paragraph here is about the method, never about the data, so it is as
# correct for someone else's study as it is for the one that produced it. Blocks
# are emitted only when the matching analysis actually ran, so nobody reads
# about a section their report does not contain.

#' A reader's guide to the analyses in a report
#'
#' Returns markdown explaining each analysis that ran: what question it asks,
#' what its output columns mean, and what the result does and does not license.
#' Only the analyses present in `res` are described.
#'
#' @param res The results list written by the workflow's analysis stage.
#'
#' @return A character vector of markdown lines.
#' @export
ap_report_guide <- function(res) {
  out <- character(0)
  add <- function(...) out <<- c(out, ..., "")

  add("# How to read this report", "",
      "This section explains what each analysis does and how to read its output. It",
      "describes only the analyses that ran for this dataset. Nothing here depends on",
      "your results, so it says the same thing for every study.")

  # --- design -----------------------------------------------------------------
  add("## Study design", "",
      "Before any test was chosen, every metadata column was classified by role:",
      "`constant`, `binary`, `categorical`, `numeric`, `high-cardinality` or",
      "`identifier`. Only `binary` and `categorical` columns can group samples.",
      "`smallest_level_n` is the size of the smallest group, and it limits what any",
      "test can detect far more than the total sample size does.", "",
      "This stage does not choose the test. It reports what the design is so the",
      "choice is made with the facts visible.")

  if (!is.null(res$scan$batch_candidates) && length(res$scan$batch_candidates) > 0L) {
    add("**Batch and technical proxies.** Columns whose names suggest a sequencing run,",
        "plate, extraction round or collection site were flagged. Test these alongside",
        "your variable of interest. If a technical variable explains variation",
        "comparable to the biological one, that belongs in the results, not in a",
        "limitations paragraph.")
  }

  if (!is.null(res$scan$repeated_measures) && length(res$scan$repeated_measures) > 0L) {
    add("**Repeated measures.** One or more columns take the same value in more than one",
        "sample. If those are subject identifiers then the samples are not independent,",
        "and any test that assumes independence will report a confidence interval that",
        "is too narrow and a p-value that is too small. Tests must use the subject as a",
        "random effect.")
  }

  if (!is.null(res$scan$confounders)) {
    add("**Candidate confounders.** Other variables were screened for association with",
        "your variable of interest, using Kruskal-Wallis for numeric columns and a",
        "chi-squared test of independence for categorical ones. A flagged variable is",
        "one to think about and report. The screen says two variables move together. It",
        "does not say which one matters, and it does not say which adjustment is right.")
  }

  # --- depth ------------------------------------------------------------------
  if (!is.null(res$depth) && !is.na(res$depth)) {
    add("## Sequencing depth", "",
        "Samples were subsampled to a common depth before the diversity metrics that",
        "need one. Rarefying discards data, and it is done because richness rises with",
        "depth, so comparing unequal libraries measures sequencing effort as much as",
        "biology. Samples below the depth were dropped, and the count of those is",
        "reported. A depth that keeps every sample is not automatically the right one:",
        "it may be so low that it discards most of the data.")
  }

  # --- alpha ------------------------------------------------------------------
  if (!is.null(res$alpha_test)) {
    add("## Alpha diversity", "",
        "Diversity within each sample, reported as Hill numbers. The naming trips people",
        "up, so read this before looking for a metric you cannot find:", "",
        "- `q0` is observed richness, the number of features present.",
        "- `q1` is the exponential of Shannon entropy. **This is Shannon diversity.**",
        "- `q2` is inverse Simpson. **This is Simpson diversity.**",
        "- `evenness` is Pielou's evenness, how equally abundance is spread.",
        "- `faith_pd` is Faith's phylogenetic diversity, the total branch length spanned.", "",
        "Hill numbers are used because they are all in units of effective number of",
        "features, so q0, q1 and q2 are directly comparable to each other. Raw Shannon",
        "entropy is in bits and is not. The three differ in how much they weight rare",
        "features: q0 counts them equally with abundant ones, q2 nearly ignores them.",
        "Disagreement between q0 and q2 is informative rather than a problem.", "",
        "**The test was chosen from the design, not picked in advance.** Two groups give",
        "a t-test with Hedges' g, or a Wilcoxon rank-sum test with Cliff's delta when the",
        "normality or equal-variance check fails. More than two groups give a one-way",
        "ANOVA with eta squared, or Kruskal-Wallis with an eta squared based on H.",
        "Covariates give a linear model and partial eta squared. Repeated measures give a",
        "linear mixed model with a random intercept per subject. The test used is named",
        "in the `test` column of every row.", "",
        "`p` is the raw p-value and `q` is that value after Benjamini-Hochberg correction",
        "across the metrics tested. Read `q`, not `p`. Read the effect size before",
        "either: at a few hundred samples a small p-value is cheap, and the effect size",
        "with its interval is what a claim rests on.")
  }

  # --- beta -------------------------------------------------------------------
  if (!is.null(res$permanova)) {
    add("## Beta diversity", "",
        "Differences in community composition between samples. Each distance answers a",
        "different question, which is why more than one is reported:", "",
        "- **Bray-Curtis** uses abundances and ignores phylogeny. Dominated by common features.",
        "- **Jaccard** uses presence and absence only, so rare features count as much as abundant ones.",
        "- **Unweighted UniFrac** uses presence and absence plus the tree, so it is sensitive to rare lineages.",
        "- **Weighted UniFrac** uses abundance plus the tree, so it is dominated by abundant lineages.", "",
        "Two metrics disagreeing is a result, not an error. It usually says whether the",
        "difference lives in which taxa are present or in how much of them there is.", "",
        "**PERMANOVA** (`vegan::adonis2`) asks whether the average composition differs",
        "between groups. `R2` is the fraction of total variation the term explains,",
        "`pseudo_F` is the test statistic, and `p` comes from permuting group labels. The",
        "default here fits every term after all the others, so a term's result does not",
        "depend on the order the terms were listed in.", "",
        "**An R2 of a few percent is normal in microbiome data and is not evidence of a",
        "weak or absent effect.** Between-sample variation is enormous. Report the R2",
        "rather than describing the effect with a word.")

    if (!is.null(res$permanova$dispersion)) {
      add("**The dispersion test is not optional here.** PERMANOVA is sensitive to unequal",
          "spread as well as to a shift in centre, so a significant result alone cannot",
          "tell you which of the two you have. `vegan::betadisper` with a permutation test",
          "asks whether the groups differ in how spread out they are. When both tests are",
          "significant, the honest report gives both numbers and does not call the result",
          "a shift in composition. When dispersion alone is significant, the groups sit in",
          "the same place and one is more variable, which is a finding in its own right.",
          "A continuous term cannot be tested this way, and is marked as such.")
    }
  }

  # --- ordination -------------------------------------------------------------
  if (!is.null(res$ordination_diagnostics)) {
    add("## Ordination plots", "",
        "The ordination panels compress the full distance matrix into two dimensions so",
        "it can be looked at. **An ordination is a picture, not a test.** The statistical",
        "answer is the PERMANOVA above. Points that look separated may not be, and points",
        "that overlap may still differ.", "",
        "The diagnostics table says whether the picture can be read as a map:", "",
        "- **PCoA** preserves distances directly. `negative_eigenvalue_fraction` is the",
        "  share of eigenvalue mass that is negative, which happens when the distance is",
        "  not Euclidean. A large value means the projection distorts the real distances.",
        "  `axes_1_2_explained` is the variation carried by the two axes drawn. When it is",
        "  low, most of the structure is in dimensions not shown.",
        "- **NMDS** preserves rank order instead, and `stress` measures how badly it had",
        "  to compromise. Below 0.1 the arrangement represents the distances well. Between",
        "  0.1 and 0.2 it is usable for broad pattern but not fine distinctions. Above 0.2",
        "  it should not be read as a map at all.", "",
        "A verdict other than `readable` or `good` means the figure needs a caveat in any",
        "caption that uses it.")
  }

  # --- composition ------------------------------------------------------------
  add("## Composition", "",
      "The taxonomic bar charts and heatmap are descriptive. They show which taxa are",
      "present and in what proportion, at the ranks given. No test is attached to them,",
      "and differences visible in a stacked bar are not evidence of anything on their",
      "own. Relative abundances sum to one, so a taxon appearing to fall may simply",
      "reflect another rising. Use the differential abundance results for claims about",
      "individual taxa.")

  # --- differential abundance -------------------------------------------------
  if (!is.null(res$da)) {
    add("## Differential abundance", "",
        "Which individual features differ between groups. Up to four methods are run on",
        "one common, pre-filtered feature set: ANCOM-BC2, ALDEx2, LinDA and MaAsLin2.", "",
        "**More than one method is run on purpose.** These methods disagree substantially",
        "on real data because they make different assumptions about compositionality,",
        "zeros and normalization. A result found by one method alone is a weaker claim",
        "than one found by several, and reporting a single method's list without saying",
        "so overstates what is known.", "",
        "Features present in too few samples are removed before testing, once, for all",
        "methods. The comparison is made against a stated reference level, which is",
        "recorded so the sign of an effect is unambiguous. Read the adjusted p-value, not",
        "the raw one.")
  }

  if (!is.null(res$concordance)) {
    add("**Method concordance.** The agreement matrix is the Jaccard overlap between each",
        "pair of methods on the features they called significant. The consensus set holds",
        "features called by at least the stated number of methods that also agree on the",
        "direction of the effect. Features called by several methods with opposite signs",
        "are flagged: those are not findings, they are a warning that the feature behaves",
        "differently under different assumptions.")
  }

  # --- screen -----------------------------------------------------------------
  if (!is.null(res$screen)) {
    add("## Exploratory screen", "",
        "Every usable metadata variable tested against alpha and beta diversity, to",
        "suggest what to look at next.", "",
        "**These are candidates, not findings.** Reporting a screen hit as a result is",
        "the garden of forking paths: with enough variables something always reaches",
        "significance. A hit here is a hypothesis to test properly, ideally in other data.", "",
        "Rows are ranked by variance explained **adjusted for degrees of freedom**, never",
        "by p-value. Raw R2 rises with the number of levels a variable has, so an",
        "unadjusted ranking puts many-level nuisance variables above real two-group",
        "effects. One Benjamini-Hochberg correction is applied across every test run, not",
        "per family. `stability` is how often a variable held its rank when samples were",
        "resampled. A hit that survives most resamples is a different claim from one that",
        "survives few, and that frequency is part of the result.")
  }

  # --- explains ---------------------------------------------------------------
  if (!is.null(res$explains)) {
    add("## Confirmatory model", "",
        "The screen tests variables one at a time. This fits them together in one model,",
        "so they compete. Each term's variance is reported after all the others, which is",
        "why a variable can look strong in the screen and weak here: the screen credited",
        "it with variation another variable also explains.", "",
        "`shared_R2` is the variance the model explains that no single term can claim.",
        "When it is large the terms are entangled and attributing the effect to one of",
        "them is not supportable from these data.")
  }

  # --- normalization ----------------------------------------------------------
  if (!is.null(res$normalization)) {
    add("## Normalization sensitivity", "",
        "Sequencing gives relative, not absolute, abundances, and there is no consensus",
        "on how to handle that. The same comparison is therefore repeated under several",
        "normalizations: total-sum scaling, cumulative-sum scaling, a centred log-ratio",
        "transform with Aitchison distance, and rarefaction.", "",
        "**The point is not to pick the best one. It is to find out whether your",
        "conclusion depends on the choice.** A result that holds under all of them is a",
        "much stronger claim than one that appears under a single normalization. When",
        "they disagree, the disagreement is the finding, and reporting one normalization",
        "without mentioning the others is not defensible.")
  }

  # --- limits -----------------------------------------------------------------
  add("## What this report does not do", "",
      "- It does not establish causation. Every test here is associational.",
      "- It does not adjust for confounders on its own. Candidate confounders are",
      "  reported; deciding what to adjust for is yours, and it is a judgement about the",
      "  study, not a statistical one.",
      "- It does not choose your hypothesis. The screen is exploratory by construction.",
      "- It does not correct a design problem. If the smallest group is tiny, or the",
      "  variable of interest tracks the sequencing batch, no analysis here repairs that.",
      "- Thresholds used anywhere in this report are choices recorded with the result,",
      "  not community standards. They are stated so you can disagree with them.")

  out
}
