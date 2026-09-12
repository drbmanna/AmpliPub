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
| `--sra-fallback` | Use `fasterq-dump` for runs ENA does not mirror, and for runs whose ENA files fail size or MD5 verification |
| `--no-metadata` | Skip the NCBI BioSample attribute lookup |

Rerunning the same command resumes. Files that already pass their checks are skipped,
and partial files continue from where they stopped.

**ENA first, SRA only as a fallback.** ENA serves the submitted FASTQ files and
publishes a size and MD5 for each, so every byte is checked against the archive.
`fasterq-dump` rebuilds FASTQ from the SRA archive, so there is no archive checksum
to compare against and the read headers are rewritten. The fallback therefore checks
the read count against ENA's `read_count` and records the run as `sra` in
`run_sources.tsv`. When a file fails verification the whole run is refetched, not the
single broken mate, so both mates come from one source and stay in the same order.

### Outputs

| File | Contents |
|---|---|
| `fastq/` | The read files, named as on ENA |
| `manifest.tsv` | QIIME 2 `PairedEndFastqManifestPhred33V2` (or single-end), one row per run |
| `run_to_sample.tsv` | Run to BioSample map, usable as QIIME 2 metadata |
| `sample_metadata.tsv` | One row per BioSample with all submitter attributes |
| `runs.tsv`, `ena_filereport.tsv` | The selected runs, and ENA's raw answer |
| `checksums.md5` | Verify with `md5sum -c checksums.md5` from the output directory |
| `run_sources.tsv` | Per run: `ena` or `sra`, and what verified it |
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

## 02_quality.py

Proposes the DADA2 truncation lengths from read quality, and stops if the truncated
reads could not merge. Runs in the QIIME 2 environment. Standard library Python 3.8 or
later.

```bash
python workflow/02_quality.py -i ~/research/baxter2016/q2/primers/trimmed.qza -o ~/research/baxter2016/q2/quality
```

DADA2 in QIIME 2 has no default truncation length. This script applies a rule fixed
before the data were seen: each read is cut just before the first position whose median
quality falls below Q30, R1 and R2 separately. It stops at the first dip even if quality
recovers later.

The truncated reads must still overlap. The floor is amplicon length + DADA2 minimum
overlap + a margin for length variation: 253 + 12 + 20 = 285 bp for V4. Below it, merging
fails and DADA2 returns a nearly empty table without an error, so the script stops and
reports both lengths.

| Option | Meaning |
|---|---|
| `--min-q Q` | Median quality a position must reach, default 30 |
| `--amplicon-len`, `--min-overlap`, `--margin` | The overlap floor, defaults 253, 12 and 20 |
| `--n N` | Reads sampled for the quality profile, default 10000 |
| `--timeout S` | Seconds before `demux summarize` is killed, default 3600 |

`qiime demux summarize` samples reads at random and takes no seed, so a rerun can shift
a median slightly. The exact table behind the decision is kept.

### Outputs

| File | Contents |
|---|---|
| `trunc_len.tsv` | `trunc_len_f`, `trunc_len_r`, the rule's settings, and the expected overlap. Read by `03_dada2` |
| `quality_profile.tsv` | Median quality and read count at every position, R1 and R2 |
| `quality.qzv` | The QIIME 2 quality plots the numbers came from |
| `quality_log.txt` | Command, rule, lengths, and the share of reads long enough to keep |

### Guards

| Guard | Why |
|---|---|
| Truncated reads must overlap by the floor | Merging otherwise collapses without an error |
| Median below the threshold at position 1 is fatal | Nothing usable would be left |
| Both quality tables must exist, positions 1..N in order, every row complete | A single-end input or a malformed table would give a wrong length silently |

## 03_dada2.py

Denoises the trimmed paired reads with DADA2 and judges the run against criteria written
down before it ran. Runs in the QIIME 2 environment. Standard library Python 3.8 or
later.

```bash
python workflow/03_dada2.py -i ~/research/baxter2016/q2/primers/trimmed.qza -t ~/research/baxter2016/q2/quality/trunc_len.tsv -o ~/research/baxter2016/q2/dada2 --threads 12
```

**There is no published pass/fail threshold for a DADA2 run.** Six sources were checked
on 2026-09-12: the DADA2 tutorial, the mothur MiSeq SOP, Kozich et al. 2013, Callahan et
al. 2016, Estaki et al. 2020, and the QIIME 2 denoising tutorial. Every one reports what
it observed and leaves the judgement to the reader. So this script hard-fails only where
the output is wrong rather than merely poor, and flags everything else with the source
quoted in the message. The criteria are written to `criteria.tsv` next to the results, so
the run and the standard it was held to travel together.

The two flag rules are the DADA2 tutorial's own words, made countable: "Outside of
filtering, there should no step in which a majority of reads are lost." (sic) and "If
most of your reads were removed as chimeric, upstream processing may need to be
revisited." Loss inside filtering is never flagged, because the tutorial excludes it.

`--p-n-threads` is always passed. QIIME 2's default is 1, so an unset value quietly runs
single-threaded. `--p-trunc-len-f`, `--p-trunc-len-r` and `--p-min-overlap` come from
`trunc_len.tsv`, so the overlap arithmetic and the run use the same numbers.

**One semantic trap**, read out of `run_dada.R` in this build rather than assumed: in the
paired track table `denoised` is `denoisedF`, the forward reads only. The drop from
`denoised` to `merged` therefore carries reverse denoising as well as merging. The
reports say so rather than calling it a merge rate.

| Option | Meaning |
|---|---|
| `-t`, `--trunc-len FILE` | `trunc_len.tsv` from `02_quality.py` |
| `--threads N` | DADA2 threads, default 1, always passed to QIIME 2 explicitly |
| `--majority F` | Fraction that counts as the tutorial's "majority", default 0.5 |
| `--allow-zero-read-samples` | Downgrade the zero-read library check to a flag. Recorded in `criteria.tsv` |
| `--timeout S` | Seconds before `denoise-paired` is killed, default 86400 |

### Outputs

| File | Contents |
|---|---|
| `table.qza`, `rep_seqs.qza` | The feature table and the ASV sequences |
| `denoising_stats.qza` | The QIIME 2 stats artifact as produced |
| `dada2_stats.tsv` | Per sample: every count, the loss at each step, and the share of input retained |
| `dada2_flags.tsv` | One row per sample and step that crossed the flag, with the counts behind it |
| `criteria.tsv` | The criteria applied, with the source for each |
| `dada2_log.txt` | Command, versions, parameters, pooled losses and every flag |

### Guards

| Guard | Why |
|---|---|
| Overlap is rechecked before the run starts | This is where the compute is spent; a bad length would otherwise cost hours and return a near-empty table |
| A library with zero reads is fatal unless overridden | Silently carrying a dead sample through the analysis is worse than stopping |
| An empty feature table is fatal | The most common catastrophic 16S failure, and it does not raise an error on its own |
| The stats table must have the paired-end columns, in order, all numeric | A single-end input or a changed format would give wrong percentages silently |
| `denoise-paired` must exit zero and write all three artifacts | An exit code alone does not prove the outputs exist |

## 04_mock.py

Compares the ASVs recovered from mock community samples with a reference of known
composition. Runs in the QIIME 2 environment. Standard library Python 3.8 or later.

```bash
python workflow/04_mock.py -b ~/research/baxter2016/q2/dada2/table.qza -r ~/research/baxter2016/q2/dada2/rep_seqs.qza -m ~/research/baxter2016/ref/HMP_MOCK.v35.fasta -s mock1,mock2,mock5,mock6,mock7 -o ~/research/baxter2016/q2/mock
```

**This stage reports. It does not pass or fail.** No published source sets an acceptance
threshold for a mock community, so the numbers are printed with published figures beside
them and the reader judges. Kozich et al. 2013 and Callahan et al. 2016 are quoted in the
log, labelled as comparisons and not as thresholds, with a note that OTUs and ASVs are
not the same unit.

The reference is not assumed to be pre-trimmed. The script finds the forward primer and
the reverse complement of the reverse primer in each record, takes what lies between
them, and collapses the results, because different strains can share one V4 sequence and
would otherwise be counted as several targets that could never all be recovered.
`--max-primer-mismatch` rescues a record whose primer site is not exact; those records are
named in the log, since a mismatched primer site can mean the template amplifies poorly
and a missing target is then not the pipeline's fault.

**Two limits, stated rather than buried.** The mismatch rate is counted against the
nearest reference of the *same length*. It is not mothur's `seq.error`, which aligns
first, so do not report it under that name. An ASV with an insertion or deletion has no
same-length reference and is counted separately instead of being quietly dropped.

| Option | Meaning |
|---|---|
| `-s`, `--mock-samples` | Comma-separated sample ids, or a file with one per line |
| `--primer-f`, `--primer-r` | Primers, default 515F and 806R |
| `--max-primer-mismatch N` | Mismatches allowed per primer site in the reference, default 1. Exact sites always win |
| `--low-depth-note N` | Log a note for any mock sample below N reads. Default 0, off, so no threshold is implied |
| `--timeout S` | Seconds before an export step is killed, default 3600 |

### Outputs

| File | Contents |
|---|---|
| `targets.tsv` | Every distinct target sequence, the reference names that produce it, and any record with no region found |
| `mock_summary.tsv` | Per sample: reads, ASVs, targets recovered exactly, spurious ASVs and their read share, mismatch rate |
| `mock_missing_targets.tsv` | Which targets were not recovered exactly, by sample |
| `mock_other_asvs.tsv` | Every ASV that is not an exact target, with its reads and distance to the nearest same-length reference |
| `mock_log.txt` | Command, versions, the reference breakdown, per-sample numbers, and the published comparisons |
