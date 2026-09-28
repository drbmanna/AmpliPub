# AmpliPub 0.0.2

## Differential abundance

* The prevalence filter is still pooled over the compared groups, on purpose: a per-group
  filter uses the group labels and can inflate the false positive rate. It now also reports
  what that costs. `ap_da()` returns a `filter_removed` table listing features dropped by the
  pooled filter that reach `prv_cut` within at least one group, with pooled and per-group
  prevalence, ranked by their highest per-group prevalence. The print method, the reader's
  guide, the methods section, the preprocessing summary and the figure legends all note these
  features when present. On Baxter this surfaces the cancer-associated *Peptostreptococcus*
  and *Porphyromonas* ASVs that the pooled cut had removed silently. Calls are unchanged.
* The methods section and the preprocessing summary now state that ANCOM-BC2 detects
  structural zeros and calls them on direction alone, without an effect or a p-value.

## Workflow

* `quality` now cross-checks `amplicon_len` against the region the configured primers cut
  from the mock reference, before the denoising run rather than after it. The primer pair
  and `amplicon_len` are independent settings that have to describe one region, and nothing
  caught them disagreeing: the overlap floor is built from `amplicon_len`, so a value left
  at V4's 253 bp while the primers amplify V3-V4 passes the overlap check, DADA2 runs to
  completion, and almost nothing merges. Groundwork for longer amplicon regions.
* The check is skipped, and says so in the log, when no mock reference is configured. An
  unverified number that looks verified is worse than one known to be unchecked.
* New `quality.len_tolerance` (default 50 bp), a chosen value and not a published
  threshold: within one region references vary by tens of bases, between regions by
  hundreds.

## Reporting fixes

* The figure legend no longer gives denoising as the reason Chao1 and ACE are absent unless
  the table really has no singletons. AmpliPub accepts any abundance table, and on an OTU
  table that kept its singletons that reason was false. `ap_alpha()` now records
  `has_singletons` from the input table, before rarefaction.
* Taxonomic composition legends name the grouping variable through its configured
  publication label rather than the raw column name.
* `04_mock.py` logged a literal `0.01%%` for the Kozich et al. error rate, because the
  logging call passes no `%`-arguments and so does no `%`-substitution.

# AmpliPub 0.0.1

First public release.

## Workflow

* Snakemake workflow from an SRA accession to a report: download with checksums, read QC,
  primer removal, truncation with an overlap check, DADA2, mock community recovery,
  taxonomy, pooling of runs per sample, filtering, tree, and diversity.
* Four conda environments built from committed lockfiles, and a Dockerfile that recreates
  them. Every run records the code commit, the environment exports and a lock check.

## R package

* Import from QIIME 2 artifacts, BIOM or TSV, with a canonical feature order.
* Metadata scan for repeated measures, batch variables and confounders, and a guard that
  refuses a grouping variable no test can use.
* Alpha diversity as Hill numbers, tested with effect sizes and bootstrap intervals.
* Beta diversity on four metrics, PERMANOVA with a dispersion test every time, pairwise
  PERMANOVA, and ordination diagnostics.
* Normalization sensitivity across TSS, CSS, CLR and rarefying.
* Differential abundance with ALDEx2, ANCOM-BC2, LinDA and MaAsLin2, and a concordance
  report.
* A screen that ranks candidate variables within each test family.
* Publication figures as PDF, SVG and 600 dpi PNG with legend files, a reader's guide, and
  a methods section generated from what ran.

## Validation

* See `VALIDATION.md`. Tested on one dataset (Baxter et al. 2016, 16S V4).

## Known issues

* The differential abundance prevalence filter is computed over the compared groups
  together and silently removes features common in one group and rare in the other. On
  Baxter it removed *Peptostreptococcus* and *Porphyromonas* ASVs associated with cancer.
  Planned for 0.0.2: list such features in the output and the report.
