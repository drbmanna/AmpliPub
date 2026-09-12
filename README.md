# AmpliPub

From raw amplicon reads to publication-grade statistics and figures, in one reproducible
project.

**Status: early development.** The first stages run today. Most of the pipeline and the
whole R statistics layer are still being built, and the interface will change.

## Why

A 16S analysis usually means QIIME 2 in one place, R in another, statistics in a third,
and plotting in a fourth. To make life simple, AmpliPub runs all of it as one project:

- Downloads reads from a single accession and checks that every file arrived intact.
- Runs read QC and QIIME 2 processing, and stops when a step fails silently.
- Runs the statistics with the checks each test needs, such as a dispersion test with
  every PERMANOVA.
- Produces figures ready to submit, and a methods draft from the analysis that ran.
- Logs every step, so the analysis can be rerun.

## How it fits together

```
accession -> FASTQ -> read QC -> QIIME 2 processing -> feature table -> statistics -> figures
```

| Stage | What runs | Status |
|---|---|---|
| Download | `workflow/00_fetch_sra.py`: one accession in, verified FASTQ, sample metadata, and a QIIME 2 manifest out | Built, tested |
| Environments | `workflow/setup_envs.sh`: the QIIME 2 amplicon 2025.7 release plus a pinned QC environment | Built |
| Read QC | `workflow/00_qc_raw.py`: FastQC and MultiQC, before and after primer removal | Built, tested |
| Primer removal | `workflow/01_primers.py`: QIIME 2 (cutadapt), anchored primers, reads counted in and out | Built, tested |
| Quality and truncation | `workflow/02_quality.py`: QIIME 2 quality profiles, with an overlap check before reads are merged | Built, tested |
| Denoising | `workflow/03_dada2.py`: QIIME 2 (DADA2), with read retention checked per library against pre-registered criteria | Built, not yet run on a full dataset |
| Mock community check | `workflow/04_mock.py`: recovered sequences compared with a community of known composition, reported against published figures | Built, not yet run on a full dataset |
| Taxonomy | `workflow/05_taxonomy.py`: QIIME 2 classifiers, with a disagreement rate between two references reported per rank | Built, not yet run on a full dataset |
| Tree | QIIME 2 | Planned |
| Statistics and figures | AmpliPub R package | Planned |

QIIME 2 does the sequence processing, and AmpliPub does not reimplement it. What AmpliPub
adds around it is getting data in from a single accession, QC reports, and checks that stop
the run on failures that otherwise pass silently: reads that no longer overlap enough to
merge, libraries that lose most of their reads, a mock community that does not come back.

Already have a table from another pipeline? Start at the feature table. The R package
takes any OTU, ASV, or gene count table plus sample metadata, whether it came from
QIIME 2, DADA2, mothur, or anything else.

The pipeline is being built one stage at a time on a public dataset (Baxter et al. 2016)
that includes mock community samples, so each stage can be checked against a known answer.
`R CMD check` runs on every push.

## Planned statistics layer

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

## Installation

The R package has nothing to install yet.

The workflow runs on Linux or WSL and needs conda:

```bash
bash workflow/setup_envs.sh
```

See [workflow/README.md](workflow/README.md) for each stage.

## License

MIT. See [LICENSE.md](LICENSE.md).
