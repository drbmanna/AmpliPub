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
