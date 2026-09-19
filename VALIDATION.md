# Validation

What has been checked in AmpliPub, against what, and with what result. Every number below
comes from a recorded run. Checks that have not been done are listed at the end, because
where a tool has not been tested is part of what it can claim.

The development dataset is Baxter et al. 2016 (Genome Medicine), SRA project PRJNA290926:
490 stool samples and 5 mock community samples, 16S V4, sequenced on MiSeq. It was run end
to end through the workflow in this repository. After filtering and rarefaction to 10,000
reads, 487 samples remain.

## 1. Agreement with QIIME 2

QIIME 2 amplicon 2025.7, whose diversity metrics come from scikit-bio, is used as an
independent reference. AmpliPub recomputes every metric in R from the same rarefied table,
so the random draw is not part of the comparison. Maximum absolute difference over all 487
samples (or all sample pairs):

| Quantity | Max abs difference |
|---|---|
| Observed features | 0 (exact) |
| Shannon entropy (bits) | 1.8e-15 |
| Pielou's evenness | 2.2e-16 |
| Faith's PD | 5.9e-08 |
| Bray-Curtis | 0 (exact) |
| Jaccard | 1.1e-16 |
| Unweighted UniFrac | 3.2e-07 |
| Weighted UniFrac | 1.9e-07 |
| PCoA proportion explained | 1e-06; axis correlations above 0.9999 |

The residuals trace to the float32 branch lengths in the Newick tree. Two conventions were
settled by this comparison rather than assumed: QIIME 2's weighted UniFrac is the raw form
(the normalised form differs by 0.38), and QIIME 2 reports Shannon entropy in bits where
vegan reports nats.

These tests are in `tests/testthat/test-crosscheck-qiime.R`. They need the QIIME 2 output,
which is research data and is not in this repository, so they skip unless the environment
variable `AMPLIPUB_QIIME_DIR` points to it.

## 2. Hand-worked reference values

Agreement on real data can hide a disagreement that only shows on some inputs. Small cases
with known answers are checked as well.

- **Unweighted UniFrac.** Eight sample pairs on two small trees, with reference values from
  scikit-bio 0.6.2. The UniFrac in `mia` (through rbiom 2.2.1) gave a different value on
  all four pairs that do not span the root (for example 0.5 where scikit-bio and phyloseq
  give 1/3). On the Baxter data nearly every pair spans the root, so the QIIME 2 comparison
  above could not see this. AmpliPub therefore keeps its own unweighted UniFrac and tests it
  against all eight scikit-bio values. Weighted UniFrac, where all three agree, comes from
  `mia`.
- **Chao1 and ACE** agree with `vegan::estimateR` to within 1e-10 over 25 random count
  vectors.

## 3. Planted truth

Synthetic data where the right answer is fixed before any method runs.

- **Differential abundance.** A table with five features enriched in one group
  (`tests/testthat/helper-fixtures.R`). The tests check that the methods recover them,
  stay specific, and report the direction of the effect the right way round.
- **ANCOM-BC2 structural zeros.** A feature planted as absent from one group is run through
  the real `ANCOMBC::ancombc2` and must come back as a structural zero with the right
  direction.
- **The screen.** A null variable with many levels must not outrank a real two-level
  effect, which is what ranking on unadjusted R2 gets wrong.

## 4. Mock communities

The Baxter data include five samples of a mock community, sequenced in the same runs as the
patients, with its reference sequences (`HMP_MOCK.v35.fasta`, 25 distinct V4 targets). The
workflow reports recovery for each sample without applying a pass or fail threshold:

| Run | Reads | Targets recovered exactly | Reads matching a target exactly | Reads not attributable to the reference |
|---|---|---|---|---|
| SRR2144132 | 750 | 2 of 25 | 2.9% | 96.7% |
| SRR2144133 | 65,372 | 22 of 25 | 94.8% | 5.1% |
| SRR2144134 | 49,751 | 17 of 25 | 61.2% | 37.8% |
| SRR2144135 | 44,603 | 3 of 25 | 11.3% | 84.4% |
| SRR2144136 | 43,818 | 3 of 25 | 10.2% | 86.4% |

SRR2144133 behaves as a mock should: 22 of 25 targets, and no ASV lies within 3 mismatches
of a reference without matching it exactly. SRR2144132 has too few reads to judge. SRR2144135 and
SRR2144136 are mostly sequences that match nothing in the reference file. Whether those two
samples used a different community, or were contaminated, has not been established, so they
are reported here as observed and not counted for or against the pipeline.

## 5. Reproducibility

- **Same inputs, run twice.** The R stage run twice on identical upstream output gave all 22
  result tables byte-identical, and 14 of 15 figures. PNG files are not byte-identical
  across R sessions (font caching), so figures are not claimed to be.
- **Same values in a different order.** Two upstream runs gave the same values with the
  features in a different order. The R stage now imposes one order on import, and 21 of 22
  tables came out byte-identical. The 22nd compares against QIIME 2's own rarefied table,
  which QIIME 2 draws without a fixed seed, so it is expected to differ.
- **Locked environments.** All four conda environments were rebuilt from the committed
  lockfiles and matched them package for package. That rebuild used a local package cache,
  so it shows the lockfiles install, not that every download link still works.
- **Provenance.** Each run records the commit of the code that ran and how many package
  files were uncommitted at the time.

## 6. Defects these checks found

Each was fixed, and each has a test that fails if it comes back.

- **ANCOM-BC2 significance.** Every result with an adjusted p-value below 0.05 was counted
  as a call, ignoring ANCOM-BC2's own pseudocount sensitivity analysis. On Baxter, cancer
  versus normal, that gave 146 calls, of which 145 failed the sensitivity analysis. With the
  package's recommended rule applied there is 1.
- **MaAsLin2 renamed features.** MaAsLin2 passes feature names through `make.names()`,
  which changed 247 of 391 ASV identifiers. Its results silently failed to match the other
  methods. Found because the method counts did not add up.
- **UniFrac and Shannon conventions** (section 1 and 2 above).
- **A false "not rarefied" caption** on tables already at one depth.
- **The screen ranked alpha and beta tests together**, although their effect sizes measure
  different variances. The smaller beta values almost never reached the top of the list, so
  their rank stability meant nothing. They are now ranked separately.
- **Effect size clamping.** The Kruskal-Wallis eta squared from `rstatix` (1.0.0 and later)
  clamps negative values to 0, which piles null effects at zero. AmpliPub computes it from
  the documented formula instead, as a stated exception to its rule of calling established
  packages.

## 7. Not yet validated

- **One dataset.** Everything above is on one V4 dataset. Other regions (V3-V4 is next) and
  other study designs have not been run.
- **Numeric benchmark against other tools.** AmpliPub has not yet been compared
  number for number with MicrobiomeAnalystR or nf-core/ampliseq on shared statistics.
- **Conclusion-level benchmark.** Whether AmpliPub's defaults lead to fewer wrong
  conclusions than other tools, on simulations with a planted answer and on published
  datasets, is planned and not done. Its success criteria will be written down before it
  runs.
- **The two unexplained mock samples** (section 4).
