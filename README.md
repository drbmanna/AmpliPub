# AmpliPub

Publication-grade statistics and figures for amplicon sequencing data.

**Status: early development. Not yet usable. The API will change.**

## Why

Getting from a feature table to a submittable figure means writing the same few hundred
lines of `vegan` and `ggplot2` code on every project. Most people write it once, write it
differently the next time, and cannot reproduce either version a year later.

The statistics have a second problem. Several of the tests that peer review expects
depend on assumptions that are rarely checked. PERMANOVA is the clearest case: a
significant result can mean the groups differ in location, or that one group is simply
more variable. Reporting it without the accompanying dispersion test leaves that
ambiguity unresolved, and it is left unresolved in a great many published analyses.

AmpliPub runs the tests, checks the assumptions they depend on, and produces figures that
are ready to submit.

## What it does

Input is a feature table (OTU, ASV, or gene counts) plus sample metadata. Source does not
matter: QIIME 2, DADA2, mothur, or anything else that yields a table.

Planned scope:

- **Data and metadata assessment.** Sparsity, depth, and library size. Detection of
  repeated measures, batch variables, and candidate confounders, which determine whether
  later stages need mixed models.
- **Normalization with sensitivity analysis.** The choice between rarefying, TSS, CSS, and
  CLR is genuinely contested. Rather than picking one silently, AmpliPub reports whether
  your conclusion survives the choice.
- **Alpha diversity.** Hill numbers rather than ad hoc index selection, with coverage
  standardization and model-based testing.
- **Beta diversity.** Dispersion testing runs with every PERMANOVA, and the output states
  plainly whether the result is a location shift. Multiple distance metrics side by side.
- **Differential abundance.** Several methods in parallel with a concordance report,
  because they are known to disagree. Effect sizes and intervals alongside adjusted
  p-values.
- **Reporting.** Figures built with `ggplot2`, `ggsci`, and `ggpubr`. A methods draft
  populated from the analysis that actually ran. A STORMS checklist.

## Scope boundaries

AmpliPub starts at the feature table. Upstream sequence processing is out of scope:
primer removal, denoising, chimera filtering, and taxonomy assignment belong to QIIME 2,
DADA2, or an equivalent, and should be done before anything here is run. Results are only
as good as that upstream work.

## Installation

Not yet. There is nothing to install.

## License

MIT. See [LICENSE.md](LICENSE.md).
