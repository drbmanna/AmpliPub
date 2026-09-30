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

### A second dataset, two regions

Chen et al. 2020 (PRJNA643648) sequenced the same ten subjects twice, at V3-V4 with 341F-805R
on 2x300 and at V4 with 515F-806R, each run carrying a ZymoBIOMICS mock. That makes the V4 arm
a comparison for the V3-V4 arm on the same subjects and the same community.

| Arm | Samples | Reads into DADA2 | Mock run | Reads | ASVs | Targets recovered | Reads matching a target exactly | Reads not attributable |
|---|---|---|---|---|---|---|---|---|
| V3-V4, 427 bp | 29 | 4,994,246 | SRR12141642 | 60,421 | 24 | 10 of 10 | 95.47% | 4.53% |
| V4, 253 bp | 11 | 1,240,260 | SRR12141640 | 58,853 | 21 | 9 of 9 | 98.81% | 1.19% |

The target counts differ because the two primer pairs collapse the reference into different
numbers of unique regions: 10 at V3-V4 where the targets are 427 and 428 bp, 9 at V4 where they
are all 253 bp. Neither is a subset of the other, so the two recovery figures are not directly
comparable to each other. Each says that its own region reproduces a known community.

Taxonomic coverage on the same subjects, as the fraction of ASVs and of reads assigned:

| Rank | V3-V4 (1,523 ASVs) | V4 (806 ASVs) |
|---|---|---|
| family | 0.983 / 0.999 | 0.906 / 0.998 |
| genus | 0.892 / 0.989 | 0.814 / 0.973 |
| species | 0.659 / 0.787 | 0.550 / 0.697 |

**The table above is not a region comparison.** The V3-V4 arm uses a classifier trained here
from Greengenes2 2024.09 with `workflow/scripts/train_classifier.py`, because no V3-V4
classifier ships with Greengenes2. The V4 arm uses Greengenes2's own V4 classifier. The two
differ in training procedure as well as in region, so those numbers are region plus classifier.

To separate them, both arms were classified again with one classifier, full-length Greengenes2
2024.09, at the same confidence of 0.7. One reference, one training procedure, region the only
difference. The criterion was fixed before the run: no region effect would be claimed if the
ASV-weighted and read-weighted coverage disagreed in direction.

| Rank | V3-V4 | V4 | Difference |
|---|---|---|---|
| family, by ASV | 0.973 | 0.908 | +6.49 pp |
| family, by read | 0.965 | 0.974 | -0.83 pp |
| genus, by ASV | 0.855 | 0.790 | +6.46 pp |
| genus, by read | 0.933 | 0.940 | -0.71 pp |
| species, by ASV | 0.619 | 0.527 | +9.19 pp |
| species, by read | 0.759 | 0.657 | +10.16 pp |

At family and genus the two weightings point in opposite directions, so by the stated criterion
no region effect is claimed at those ranks. At species both agree and V3-V4 is ahead by 9 to 10
percentage points. **What this dataset supports is a species-level difference, nothing more.**
The larger genus gap in the first table, 0.892 against 0.814, does not survive holding the
classifier constant, so most of it was the classifier and not the region.

Read-weighted coverage is already between 0.93 and 0.97 at genus for both regions, so abundant
taxa are named from either one and the ASV-weighted gap sits in the rare tail. That is an
explanation offered for the pattern, not something measured here.

Both arms were produced by one commit with no uncommitted changes, and both reproduced when
re-run: the V4 arm returned every mock and coverage figure identical to all printed decimals,
and the V3-V4 arm moved by three ASVs out of 1,526 with coverage unchanged in the fourth
decimal place.

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
- **3' primer readthrough.** Only 5' anchored primers were trimmed. When a read is longer
  than its amplicon, sequencing runs past the far primer and reads it, and nothing removed
  that. On the V4 arm above, 97.69% of R1 and 78.13% of R2 carried the far primer, so the ASVs
  were 273 bp against a 253 bp amplicon and ended in the reverse complement of 806R, matching
  base for base including every ambiguity code. The mock reported 0 of 9 targets recovered,
  the ASV sequences were not amplicons, their lengths were wrong, and classification ran on
  sequences carrying 20 bp of primer. After the fix that arm recovers 9 of 9.
  This escaped every earlier check because validation had only ever used Baxter, where the
  reads are 251 bp and the amplicon 253, so readthrough was impossible. The condition is
  read length greater than amplicon length, which says in advance which datasets are affected,
  and the three arms now measured span it: 97.69% of R1 at V4 on 2x300, 0.15% on Baxter,
  0.11% at V3-V4 where the 427 bp amplicon cannot be crossed by a 300 bp read.
- **Two settings for one region, with one of them hardcoded.** The mock stage did not receive
  the configured primers and fell back to 515F-806R whatever the run used, so on V3-V4 the
  reference was cut to 253 bp against 427 bp ASVs and every field of the summary was a
  well-formed zero. The same defect class as `amplicon_len` disagreeing with the primer pair.
  The stage now refuses when there are ASVs, no target is recovered, and no ASV shares a length
  with any target, which is a region mismatch rather than a property of the library. All three
  conditions are required, because zero recovery on its own can be a genuinely poor mock and
  must still be reported as one.
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

- **Two datasets, two regions, one study design.** The statistics and figures are validated
  on Baxter alone. The workflow has now been run on a second dataset at two regions
  (section 4), but that dataset is ten healthy volunteers with no group contrast, so its R
  stage was deliberately not run and it tests the upstream workflow only. No other study
  design has been run.
- **Region effect above species level.** With one classifier on both arms, V3-V4 and V4 differ
  at species by 9 to 10 percentage points, but at family and genus the ASV-weighted and
  read-weighted figures disagree in direction, so no claim is made at those ranks (section 4).
- **Why region-matched classifiers beat full-length is not established.** That they do is
  measured, on three arms and against a pre-registered rule. The reason is not. The obvious
  explanation, that a full-length model classifying a short fragment has flatter posteriors and
  so fails the confidence threshold more often, predicts a smaller penalty on a longer amplicon
  and the opposite was found: the genus penalty is 3.68 pp by ASV on the 427 bp V3-V4 amplicon
  against 2.36 pp on the 253 bp V4 amplicon of the same subjects. The penalty also varies by
  dataset at fixed amplicon length, 9.40 pp on Baxter against 2.36 pp here, both at 253 bp.
  AmpliPub therefore keeps region-matched classifiers on the evidence and offers no mechanism.
- **Readthrough trimming is not exhaustive.** `primers.readthrough_overlap` defaults to 10, so
  at least 10 bases of the far primer must match. A read that sequenced only a few bases into
  that primer keeps them. On the V4 arm about 80 of 806 ASVs remain in a 268 to 293 bp band
  above the 253 bp amplicon for this reason. The trim removes the common case, not every case.
- **The retention floor is judged against an unseeded estimate.** `quality.min_reach` is
  checked against a quality profile that `qiime demux summarize` draws without a seed, so the
  retained fraction shifts slightly between runs: 0.8087 and 0.8129 on two runs of the same
  V4 data. Far from the 0.75 floor this changes nothing, and it changed nothing here. Near the
  floor it could change whether the floor is met, and so change the truncation position and
  everything downstream. The reported number is disclosed in `trunc_len.tsv`; the
  non-determinism of the verdict is not yet addressed.
- **Numeric benchmark against other tools.** AmpliPub has not yet been compared
  number for number with MicrobiomeAnalystR or nf-core/ampliseq on shared statistics.
- **Conclusion-level benchmark.** Whether AmpliPub's defaults lead to fewer wrong
  conclusions than other tools, on simulations with a planted answer and on published
  datasets, is planned and not done. Its success criteria will be written down before it
  runs.
- **The two unexplained mock samples** (section 4).
- **Features removed by the differential abundance prevalence filter.** The 10% prevalence
  filter is computed over the compared groups together. On Baxter, cancer versus normal, it
  removed *Peptostreptococcus* and three *Porphyromonas* ASVs present in 17-21% of cancer
  and 1-3% of normal samples, taxa the original study associates with cancer.
  *Parvimonas micra* (26% versus 5%) passed the filter and is called by three of
  four methods. Since 0.0.2 these features are reported rather than dropped in silence:
  `ap_da()` returns a `filter_removed` table giving pooled and per-group prevalence for
  every feature that reaches the cut in at least one group, warns when the table is not
  empty, and carries the count into the report, the methods section and the figure legends.
  On Baxter it lists 78 features, including the four above.
  The filter itself is unchanged and still pooled, on purpose: a per-group filter uses the
  group labels, so it is not independent of the test under the null and can inflate the
  false positive rate. Disclosure is not a fix. The removed features are still not tested,
  and deciding what to do with them is left to the reader.
