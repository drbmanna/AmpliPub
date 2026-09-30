# AmpliPub 0.0.3

## Workflow

* `primers` now trims the far primer from the 3' end of every read, in a second cutadapt pass
  using the reverse complement of the opposite primer. Only 5' anchored primers were trimmed
  before, so a read longer than its amplicon ran past the far primer, read it, and kept it. On
  a V4 dataset sequenced at 2x300 this affected 97.69% of R1 and 78.13% of R2: the ASVs came
  out 273 bp against a 253 bp amplicon, ending in the reverse complement of 806R, and the mock
  comparison reported 0 of 9 targets recovered. With the pass in place the same data recovers
  9 of 9. The condition is read length greater than amplicon length, so datasets where the
  read is shorter than the amplicon are unaffected: 0.15% of R1 on Baxter at 251 bp reads and
  a 253 bp amplicon, 0.11% at V3-V4 where the amplicon is 427 bp.
* New `primers.trim_readthrough` (default `true`) and `primers.readthrough_overlap`
  (default 10). The default overlap is not cutadapt's own 3, which would match by chance in
  about one read in 64 and shorten ASVs. It is a chosen value, not a published threshold, and
  it is not exhaustive: a read carrying fewer than 10 bases of the far primer keeps them.
* `quality` now caps the truncation position by how many reads survive it. The position came
  from quality alone, the first cycle whose median drops below `min_q`, and DADA2 then discards
  every shorter read without saying so. On a 2x300 V3-V4 run the quality rule chose position
  300, which only 2.2% of R1 reads reach, so the run would have completed on a fortieth of the
  forward data and reported healthy quality. The cap moves it to 284, where 95.7% survive.
* New `quality.min_reach` (default 0.75), adapted from ampliseq's `trunc_rmin`, which uses
  retention to pick a cutoff where this checks the cutoff `min_q` already picked. A chosen
  value, not a published threshold. `trunc_len.tsv` now records the quality position, the cap,
  and the retained fraction for each read, so the decision is inspectable.
* The `mock` stage now receives the configured primers. It defaulted to 515F-806R whatever
  region the run used, so on V3-V4 the reference was cut to 253 bp against 427 bp ASVs and
  every field of the summary was a well-formed zero. It also now refuses when there are ASVs,
  no target is recovered, and no ASV shares a length with any target, which is a region
  mismatch rather than a property of the library. All three conditions are required, because
  zero recovery alone can be a genuinely poor mock and must still be described as one.
* `quality` now names the case where `amplicon_len` exceeds the region the primers cut by
  exactly the combined primer length. That is the primers counted twice, not variation between
  taxa: published insert sizes usually include the primers where `amplicon_len` means the
  region without them. A configured 465 bp for 341F-805R against a measured 427 bp differs by
  38, which is exactly the two primer lengths, and a 50 bp tolerance swallowed it.
* `fetch` now records `library_name` in the run map. A BioSample is not always one library, and
  a submission registering one BioSample per study leaves `sample_accession` and `sample_title`
  identical across every run, so pooling on either merges unrelated samples in silence. The
  field is optional in ENA and is requested without being required, since ENA omits a column
  entirely when no run in the study carries it.
* New `workflow/scripts/train_classifier.py`, for regions with no pre-trained classifier.
  Greengenes2 ships one for V4 and none for V3-V4. The script extracts reads from a reference
  with the configured primers and fits a naive Bayes classifier on them, so the classifier is
  tied to the primer pair it was built for.

## Validation

* `VALIDATION.md` adds a second dataset at two regions, the readthrough defect and what it
  cost, and four limitations: one study design only, the region effect not separable from the
  classifier effect, readthrough trimming not being exhaustive, and the retention floor being
  judged against an unseeded estimate.

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
