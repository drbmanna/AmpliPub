# Workflow scripts

Staged scripts that run upstream of the AmpliPub R package, in WSL or Linux. They are
not part of the package build.

## Setup

```bash
bash workflow/setup_envs.sh
```

Creates two conda environments. `qiime2-amplicon-2025.7` comes from the official QIIME 2
release file. `amplipub-qc` holds FastQC and MultiQC, pinned in `envs/qc.yml`. The QC
tools stay out of the QIIME 2 environment so the release environment is never modified.
Rerunning is safe: existing environments are only checked, not rebuilt.

## 00_fetch_sra.py

Downloads raw FASTQ for a public accession and writes a QIIME 2 manifest. Standard
library Python 3.8 or later, no extra packages.

```bash
python workflow/00_fetch_sra.py PRJNA290926 -o ~/research/baxter2016/raw
```

The accession can be a BioProject (`PRJNA…`, `PRJEB…`), a study (`SRP…`, `ERP…`), a
BioSample (`SAMN…`) or a run (`SRR…`). Pass one accession only.

| Option | Meaning |
|---|---|
| `--dry-run` | Resolve and plan, print the size, download nothing |
| `--runs FILE` | Fetch only the run accessions listed in FILE, one per line |
| `--expect-runs N` | Fail unless the accession resolves to exactly N runs |
| `--threads N` | Parallel downloads, default 4 |
| `--retries N` | Attempts per file, default 5 |
| `--sra-fallback` | Use `fasterq-dump` for runs that ENA does not mirror |
| `--no-metadata` | Skip the NCBI BioSample attribute lookup |

Rerunning the same command resumes. Files that already pass their checks are skipped,
and partial files continue from where they stopped.

### Outputs

| File | Contents |
|---|---|
| `fastq/` | The read files, named as on ENA |
| `manifest.tsv` | QIIME 2 `PairedEndFastqManifestPhred33V2` (or single-end), one row per run |
| `run_to_sample.tsv` | Run to BioSample map, usable as QIIME 2 metadata |
| `sample_metadata.tsv` | One row per BioSample with all submitter attributes |
| `runs.tsv`, `ena_filereport.tsv` | The selected runs, and ENA's raw answer |
| `checksums.md5` | Verify with `md5sum -c checksums.md5` from the output directory |
| `fetch_log.txt` | Command, versions, queries and progress, appended on every run |

**One manifest row per SRA run, not per sample.** An SRA run (`SRR…`) is one sequenced
library, not a MiSeq flowcell run. A resequenced sample has several SRA runs, possibly
from different flowcells. Keeping them apart keeps QC and read retention visible per
library and avoids pooling reads before denoising. Sum them into samples after DADA2:

```bash
qiime feature-table group --i-table table.qza --p-axis sample --p-mode sum \
  --m-metadata-file run_to_sample.tsv --m-metadata-column sample_accession \
  --o-grouped-table table-by-sample.qza
```

### Guards

Each one stops the run with a message, and each has a test in `tests/`.

| Guard | Why |
|---|---|
| Zero runs is fatal | ENA answers unknown or private IDs with an empty table and HTTP 200 |
| Comma lists are rejected | ENA answers a list of accessions with an empty table too |
| Size and MD5 are checked before a file gets its final name | A truncated `.fastq.gz` can still decompress partway and look fine |
| A lone `_1` or `_2` is fatal | Treating a broken pair as single-end silently changes the analysis |
| Mixed paired and single-end runs are fatal | One manifest cannot hold both |
| `--expect-runs` mismatch is fatal | Catches a wrong or changed accession |
| The manifest row count must equal the run count | Checked at the output boundary |

Instrument labels are reported and **never used to filter**. SRA labels can be wrong:
467 of the 544 Baxter 2016 runs say `454 GS` but are Illumina MiSeq 2 x 251.

Runs that list an extra unpaired file keep `_1` and `_2` and log the extra file.

### Tests

```bash
python -m pytest workflow/tests -q
```

Offline, using real ENA and NCBI responses saved as fixtures. They show the guards fire.
They do not prove the downloads work against the live servers; the live run does that.

## 00_qc_raw.py

Runs FastQC on every FASTQ in a directory and combines the reports with MultiQC. Runs in
the `amplipub-qc` environment (see Setup). Standard library Python 3.8 or later.

```bash
python workflow/00_qc_raw.py -i ~/research/baxter2016/raw/fastq -o ~/research/baxter2016/qc_raw
```

| Option | Meaning |
|---|---|
| `--threads N` | Files FastQC processes at once, default 4 (about 250 MB RAM each) |
| `--timeout S` | Seconds before FastQC or MultiQC is killed, default 7200 |
| `--env NAME` | Conda environment holding FastQC and MultiQC, default `amplipub-qc` |

### Outputs

| File | Contents |
|---|---|
| `multiqc_report.html` | All FastQC reports in one page |
| `qc_summary.tsv` | One row per file: read count, each FastQC module result, and two lists, `to_check` and `expected_for_amplicons` |
| `fastqc/` | The per-file FastQC reports |
| `qc_log.txt` | Command, versions and results, appended on every run |

**Four FastQC modules fail on nearly every amplicon library, and that is expected.**
Per base sequence content, per sequence GC content, sequence duplication and
overrepresented sequences all assume random fragments of one genome. Amplicon reads all
start with the same primer and come from a mixed community. These go to
`expected_for_amplicons`. Everything else that is not PASS goes to `to_check`: quality,
N content, length and adapter content.

The script reports what FastQC found and does not stop on data quality. Decisions based
on quality, such as truncation lengths, belong to `02_quality`.

### Guards

Each one stops the run with a message, and each has a test in `tests/`.

| Guard | Why |
|---|---|
| No FASTQ files is fatal | FastQC given no files opens its GUI and waits for ever |
| An empty FASTQ is fatal | Nothing downstream should run on a zero-byte file |
| Two inputs that map to one report name are fatal | `x.fastq` and `x.fastq.gz` would overwrite each other's report |
| Every tool call has a time limit, and stdin is closed | A hung tool is killed along with anything it started |
| One report per input, naming that input | A missing or mismatched report means FastQC skipped a file |
| R1 and R2 of a run must hold the same number of reads | A pair that no longer matches breaks merging later |
| MultiQC must write its report | Checked at the output boundary |

## 01_primers.py

Removes primers with `qiime cutadapt trim-paired` and reports, per sample, how many reads
carried them. Runs in the QIIME 2 environment. Standard library Python 3.8 or later.

```bash
python workflow/01_primers.py -i ~/research/baxter2016/q2/demux.qza -o ~/research/baxter2016/q2/primers
```

Cutadapt always runs. If the reads carry no primers, nothing is cut and every read comes
out unchanged, and the report says so.

| Option | Meaning |
|---|---|
| `--primer-f`, `--primer-r` | Primers, default 515F `GTGCCAGCMGCCGCGGTAA` and 806R `GGACTACHVGGGTWTCTAAT` |
| `--discard-untrimmed` | Drop pairs where no primer was found. Only for reads that carry primers |
| `--cores N` | CPU cores for cutadapt, default 4 |
| `--timeout S` | Seconds before cutadapt is killed, default 7200 |

Two settings keep primer-free reads intact. Primers are anchored (`^`), so only a full
primer at the start of a read is cut; without the anchor, cutadapt also cuts a partial
match of 3 or more bases. Untrimmed pairs are kept by default; discarding them would drop
every read of primer-free data.

### Outputs

| File | Contents |
|---|---|
| `trimmed.qza` | Reads for the next stage |
| `primer_summary.tsv` | One row per sample: pairs in, R1 and R2 with primer, pairs out, and percentages |
| `cutadapt_report.log` | Cutadapt's full report for every sample |
| `primers_log.txt` | Command, versions, and a one-line verdict: primers absent, present, or mixed |
| `qc_trimmed/`, `trimmed_fastq/` | Only when at least one read was trimmed: the trimmed reads exported, and `00_qc_raw.py` run on them |

**QC after trimming runs only when trimming changed something.** If no read was trimmed,
the reads are identical to the raw reads, the raw-read QC from `00_qc_raw.py` applies,
and the log says the second QC was skipped. `--qc-env` names the QC environment,
default `amplipub-qc`.

### Guards

| Guard | Why |
|---|---|
| Primers must be IUPAC bases, at least 10 | A typo would silently match nothing |
| The input must be a paired demux artifact, every sample with both reads | Checked from the artifact's MANIFEST |
| One cutadapt report per sample, matched by file name to the MANIFEST | QIIME's output can run into cutadapt's report; a parser that missed those would drop samples silently |
| Pairs out must equal pairs in when untrimmed reads are kept | Any loss means something other than primer removal happened |
| The trimmed artifact must hold the same samples | Checked at the output boundary |
| A mix of trimmed and untrimmed samples is flagged | Usually a wrong primer, or libraries trimmed differently |
