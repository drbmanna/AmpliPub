# Workflow scripts

Staged scripts that run upstream of the AmpliPub R package, in WSL or Linux. They are
not part of the package build.

## Run everything

One Snakemake command runs from raw reads to the report: download, QC, import, primers,
quality, DADA2, the mock check, taxonomy, pooling runs into samples, filtering, the tree,
diversity, the AmpliPub analysis and an HTML report. Each stage below is called unchanged.

```bash
bash workflow/setup_envs.sh
cp workflow/config/config.template.yaml ~/my_study.yaml   # edit paths and checksums; keep it outside the repo
mkdir -p ~/runs/my_study && cd ~/runs/my_study
conda run -n amplipub-snakemake snakemake -s /path/to/AmpliPub/workflow/Snakefile \
  --configfile ~/my_study.yaml --sdm conda --cores 12
```

- **Dry run first.** Add `-n` to print every job without running anything.
- **Run from a working directory outside the repository.** Snakemake keeps its state in
  `.snakemake/` in the directory it is run from.
- **Every input the workflow did not produce is checked.** The metadata, the mock
  reference and the classifier are verified against the sha256 in the config before any
  stage uses them.
- **Rerunning is safe.** Finished jobs are skipped. Pointing `fetch.outdir` at an earlier
  download resumes it instead of fetching again.
- **The config is validated before anything runs**, against `workflow/config/schema.yaml`.
- **What a run leaves behind:** each stage's own outputs and log, a Snakemake log and
  benchmark per job, `amplipub/tables/` and `amplipub/figures/`, `report/report.html`,
  and `provenance/`. Provenance holds the git commit, every conda environment exported,
  `sessionInfo()`, the resolved config and the sha256 of every analysis input.

The R stages run in the `amplipub-r` environment. Every run first installs AmpliPub into
it from this checkout (rule `install_amplipub`, recorded in
`provenance/amplipub_install.tsv`), so the package that ran is the commit in
`provenance/git_commit.txt`. The run also compares each environment with its lockfile and
writes any difference to `provenance/env_lock_<env>.tsv`, with a warning in the log. A
difference does not stop the run.

## Setup

```bash
bash workflow/setup_envs.sh
```

Creates four conda environments from the lockfiles in `envs/`:
`amplipub-qiime2-2025.7`, `amplipub-qc` (FastQC and MultiQC), `amplipub-snakemake` and
`amplipub-r`. The QC tools stay out of the QIIME 2 environment so the release environment
is never modified.

A lockfile (`conda list --explicit --md5`) names every package, dependencies included,
with its exact build and md5, so no solver runs. The `*.yml` files pin only the packages
we name. On 2026-09-13 the solver filled that gap with an rbiom that mia 1.18.0 could not
load, which is why the lockfiles exist. The lockfiles are for linux-64 only.

- Rerunning is safe. An existing environment is compared with its lockfile, and the
  script stops if they differ.
- `--rebuild` removes every environment and recreates it from its lockfile.
- `--from-spec` solves the `*.yml` files instead, for a deliberate upgrade (and on other
  platforms). Test the result, then re-lock and commit:
  `conda list -n ENV --explicit --md5 > workflow/envs/ENV.lock`.

## 00_demux.py

Demultiplexes EMP-protocol paired reads and accounts for every read. Runs in the QIIME 2
environment. Standard library Python 3.8 or later.

```bash
python workflow/00_demux.py -i ~/research/run1/emp -m ~/research/run1/metadata.tsv --barcode-column barcode-sequence -o ~/research/run1/demux
```

**Already demultiplexed?** Most public data, including anything from SRA, arrives as one
FASTQ pair per sample. There is nothing to demultiplex and this stage is not the one you
want: import with a manifest, as `00_fetch_sra.py` prints at the end of its log.

**The orientation trap.** `--p-rev-comp-barcodes` and `--p-rev-comp-mapping-barcodes`
both default to False, and getting either wrong does not raise an error. It assigns
almost nothing and the run continues to a near-empty feature table that reads as a failed
experiment rather than a wrong flag. This is the most common catastrophic failure in 16S
processing. When the assigned fraction falls below `--min-assigned`, this stage tries the
other three combinations, reports what each would assign, and stops with the answer.
Nothing downstream is written.

**Read accounting is the point.** A demultiplexing step that does not tell you how many
reads it threw away is not a quality control step. Assigned plus unassigned must equal the
reads in `barcodes.fastq.gz`, and more assigned than input is fatal because the three EMP
files may not come from one run. The unassigned fraction is a measurement, not waste: it
carries barcode quality, index hopping and contamination from other libraries on the run.

The `ErrorCorrectionDetails` artifact that QIIME 2 emits and nobody reads is summarised:
records, corrections applied, and reads left with no sample.

| Option | Meaning |
|---|---|
| `-i`, `--input DIR` | Directory holding `forward.fastq.gz`, `reverse.fastq.gz`, `barcodes.fastq.gz` |
| `--barcode-column NAME` | Metadata column holding the barcodes, default `barcode-sequence` |
| `--rev-comp-barcodes`, `--rev-comp-mapping-barcodes` | The two orientation flags |
| `--no-golay` | Switch off 12nt Golay correction, which QIIME 2 has on by default |
| `--min-assigned F` | Below this the orientation is checked, default 0.5, our choice |
| `--no-orientation-check` | Do not try the other orientations on low assignment |

### Outputs

| File | Contents |
|---|---|
| `per_sample_sequences.qza` | The demultiplexed reads |
| `error_correction.qza` | The barcode error correction detail, as QIIME 2 produced it |
| `demux_counts.tsv` | Per sample: barcode, reads, and the share of the input |
| `demux_log.txt` | Command, versions, the read accounting, corrections and every warning |

### Guards

| Guard | Why |
|---|---|
| Two samples sharing a barcode is fatal | They cannot be told apart |
| Mixed barcode lengths or non-DNA values are fatal | Both mean the wrong column was named |
| Golay on for barcodes that are not 12 nt warns | Golay correction is for 12nt barcodes |
| Assigned plus unassigned must equal the input | Otherwise the accounting cannot be trusted |
| Low assignment triggers the orientation check | A near-empty table is the failure this prevents |
| A missing EMP file names the file | And says that demultiplexed data wants a manifest instead |

---

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
| `manifest.tsv` | QIIME 2 `PairedEndFastqManifestPhred33V2` (or single-end), one row per run. Paths start with `$PWD`, so the download can be moved; import it from this directory |
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

## 05_taxonomy.py

Classifies the ASVs against a reference taxonomy and reports how deep the names go. Runs
in the QIIME 2 environment. Standard library Python 3.8 or later.

```bash
python workflow/05_taxonomy.py -r ~/research/baxter2016/q2/dada2/rep_seqs.qza -b ~/research/baxter2016/q2/dada2/table.qza -c gg2=~/research/ref/classifiers/gg2-2024.09-v4.qza -c silva=~/research/ref/classifiers/silva-138-99-human-stool-weighted.qza -o ~/research/baxter2016/q2/taxonomy --n-jobs 4
```

**One reference is the default.** Pick it, pin its version, name it in the methods, and
report `taxonomy_coverage.tsv`: the share of ASVs and of reads that get a name at each
rank. That is the number a reader needs.

**A second classifier is a diagnostic, not a better default.** It tells you, once, how
much of your naming depends on the reference you chose. Two things to know before reading
the output. It measures whether the two labels are the same string, not whether the two
references place an organism differently: Greengenes2 writes GTDB names, so the phylum
SILVA calls `Firmicutes` it calls `Bacillota_A_368345`, and that counts as a different
label while naming the same clade. On the real GG2-vs-SILVA run here that put the phylum
figure at 85%, essentially all of it nomenclature. And running two references routinely
reintroduces exactly the naming problem that choosing one removes at the source.

No synonym table is applied. A partial one would silently miscount every pair of
vocabularies it did not know, and claiming agreement we cannot establish is worse than a
number that says plainly what it measures.

**The ceiling is the fragment, not the classifier.** 253 bp of V4 does not carry
species-level information for many taxa. In the mock reference used by this project,
*S. aureus.1* and *S. epidermidis.1/2* have identical V4 sequences, so no method separates
them here. A low species-level rate is a property of the amplicon and the script says so
in its log rather than leaving it to be read as a failure.

**Mock calls are listed, not scored.** With `--mock-targets` pointing at `04_mock.py`'s
`targets.tsv`, the script writes one row per reference target that is present as an exact
ASV, with each classifier's call beside it. It does not mark them right or wrong: turning
a reference name like `B.vulgatus.1` into a taxon string to compare against would mean
inventing a mapping the script has no source for. Twenty-five rows is a table a person
reads.

Classifiers are **not** shipped with this repository. Get them from the QIIME 2 Library
data resources, or build one for your own region with RESCRIPt. Verify the checksum the
Library publishes, and keep the file: the log records the artifact's UUID so a result can
be traced back to the exact classifier that produced it.

| Option | Meaning |
|---|---|
| `-c`, `--classifier NAME=PATH` | A trained classifier and the name to report it under. Repeat for more than one |
| `-b`, `--table` | Feature table, so coverage is also weighted by reads and not only by ASV |
| `--mock-targets FILE` | `targets.tsv` from `04_mock.py` |
| `--n-jobs N` | `classify-sklearn` jobs, default 1, always passed explicitly because that is QIIME 2's own default |
| `--confidence F` | `classify-sklearn` confidence, default 0.7, which is also QIIME 2's |
| `--timeout S` | Seconds before a classify step is killed, default 86400 |

### Outputs

| File | Contents |
|---|---|
| `taxonomy_<name>.qza` | The classification artifact from each classifier |
| `taxonomy_coverage.tsv` | Per classifier and rank: ASVs named, and the share of ASVs and of reads |
| `taxonomy_label_differences.tsv` | Only when two or more classifiers are given. Per pair and rank: how many features both named, how many carry a different label, and how many only one named. A string comparison, not a placement conflict |
| `taxonomy_calls.tsv` | One row per ASV, one pair of columns per classifier: the taxon and its confidence |
| `taxonomy_mock_calls.tsv` | Reference target, its exact-match ASV, and each classifier's call |
| `taxonomy_log.txt` | Command, versions, each classifier's UUID, coverage, any label differences, and the amplicon-ceiling note |

### Guards

| Guard | Why |
|---|---|
| Every classifier must really be a `TaxonomicClassifier` | The wrong artifact type would otherwise be found out hours into a run |
| Every ASV must get a taxonomy row | A silent join failure would drop features from the report rather than raising |
| Nothing placed at domain level is fatal | That is a broken classifier or the wrong region, not a hard dataset |
| Every ASV in rep-seqs must be in the table, when one is given | Guards against mixing outputs from two different runs |
| One classifier says so, instead of printing a comparison of nothing | An empty comparison table would read as agreement |

## 06_resolution.py

Measures what an amplicon region can and cannot tell apart, from a reference database and
a primer pair. No QIIME 2, no conda, no network, standard library only.

```bash
python workflow/06_resolution.py -m ~/research/baxter2016/ref/HMP_MOCK.v35.fasta -p v4=GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT --full-length -o ~/research/baxter2016/resolution
```

**Why this exists.** Every 16S study picks a region and then reports species names. Almost
none checks whether the chosen region can distinguish the species it is naming. Two
organisms with an identical sequence over the amplified region are not hard to classify,
they are impossible to classify apart, and no classifier, model or database fixes that.
This stage measures it, and it can be run before sequencing rather than after.

It reads both ways round. Choosing a region: run several primer pairs against one
reference and compare counts. Reading a result: the ambiguity sets say which species-level
calls are safe and which are one arbitrary pick from several equal candidates.

**On the reference this project uses**, 32 records of the HMP mock:

| region | records with a region | distinct sequences | collapse groups | records collapsed |
|---|---|---|---|---|
| V4, 515F/806R, 252-254 bp | 32 | 25 | 5 | 12 |
| full length, 514-542 bp | 32 | 31 | 1 | 2 |

V4 cannot separate *S. aureus* from *S. epidermidis*, nor three *B. vulgatus* variants,
nor two each of *C. beijerinckii*, *E. faecalis* and *P. aeruginosa*. At full length only
the two *E. faecalis* variants remain merged, and those are the same species. That is the
concrete version of "longer reads resolve more", with numbers instead of assertion.

**An ambiguity set is a property of the reference and the region, not of your samples.**
If a set holds three species and your sample contains one of them, the call is still
ambiguous, because nothing in the amplicon says which. Report the set.

**A region this reference cannot answer for is reported as such, not as a region that
resolves nothing.** Those are different facts and conflating them would be a wrong answer
rather than a missing one. When no record yields a region, the log says which primer site
is absent and in which orientation, because a reference trimmed to start inside the
amplicon has lost its forward site and that says nothing about the region. One failing
region never aborts a comparison; only every region failing is fatal.

| Option | Meaning |
|---|---|
| `-p`, `--primers NAME=FWD,REV` | A region and its primer pair. Repeat to compare several |
| `--full-length` | Also treat the whole record as a region, for comparison |
| `-t`, `--taxonomy FILE` | Two-column id and taxon table, to resolve by rank and name the ambiguity sets |
| `--max-primer-mismatch N` | Mismatches allowed per primer site, default 1. Exact sites always win |

### Outputs

| File | Contents |
|---|---|
| `resolution_summary.tsv` | Per region: records, distinct sequences, collapse groups, lengths, and per-rank resolved counts when a taxonomy is given |
| `resolution_groups.tsv` | Every set of records that share one sequence, no taxonomy needed to read it |
| `resolution_ambiguity_sets.tsv` | Per region and rank: the taxa that collapse together, and the records behind them |
| `resolution_no_region.tsv` | Records yielding no region, per region |
| `resolution_log.txt` | Command, versions, per-region and per-rank numbers, and the caveat |

### Guards

| Guard | Why |
|---|---|
| A record yielding no region is reported, never dropped | A mismatched primer site can mean poor amplification, a different fact from being indistinguishable |
| A failing region is reported with the missing primer named | "No region found" alone blames the region for a property of the reference |
| One failing region does not abort the comparison | Comparing regions is the entire purpose |
| Every region failing is fatal | Nothing was measured, so there is no result to report |
| Taxonomy ids that match no reference record are fatal | Otherwise every rank would silently report zero resolved |

## 07_collapse.py

Pools the sequencing runs of each sample and proves no reads were lost doing it. Runs in
the QIIME 2 environment. Standard library Python 3.8 or later.

```bash
python workflow/07_collapse.py -b ~/research/baxter2016/q2/dada2/table.qza -r ~/research/baxter2016/raw/run_to_sample.tsv --metadata ~/research/baxter2016/ref/metadata.tsv --expect-samples 495 -o ~/research/baxter2016/q2/collapsed
```

**A table from `03_dada2` has one column per run, not per sample.** A resequenced sample
is counted twice, its depth halved and its diversity measured on a fraction of its reads,
and nothing about that raises an error. Every per-sample number downstream is wrong until
the runs are pooled, which is why this stage exists and why it checks its own arithmetic
rather than trusting `--p-mode sum`.

**Sample ids are not rewritten silently.** An id with whitespace is not a usable QIIME 2
sample id, but quietly renaming somebody's samples is worse than stopping, so the default
is to fail and name the offenders. `--sanitize-ids` is an explicit opt-in that writes the
before and after map and refuses if sanitizing would collide with an existing id.

| Option | Meaning |
|---|---|
| `-r`, `--run-map FILE` | Run to sample map, e.g. `run_to_sample.tsv` from `00_fetch_sra.py` |
| `--group-column NAME` | Run-map column holding the sample id, default `sample_title` |
| `--metadata FILE` | Study metadata to attach to the pooled samples |
| `--expect-samples N` | Fail unless pooling gives exactly this many samples |
| `--sanitize-ids` | Replace whitespace in sample ids with underscores, recorded in a map |

### Outputs

| File | Contents |
|---|---|
| `table_by_sample.qza` | The pooled feature table |
| `runs_per_sample.tsv` | Per sample: how many runs, which ones, and the reads they hold |
| `sample_metadata.tsv` | The study metadata joined to the pooled samples, with `in_study_metadata` marking what matched |
| `grouping.tsv`, `sanitized_ids.tsv` | The grouping actually used, and any id that was changed |
| `collapse_log.txt` | Command, versions, the read accounting and every mismatch |

### Guards

| Guard | Why |
|---|---|
| The read total must be identical before and after | Pooling moves reads between columns; it must not create or destroy any |
| Every sample must equal the sum of its runs | The overall total can match while one sample is wrong |
| A run in the table with no map entry is fatal | `feature-table group` would drop it silently |
| Unusable sample ids stop the run | Renaming someone's samples without asking is worse than stopping |
| Metadata matching nothing is fatal | Otherwise every downstream group test runs on empty columns |

---

## 08_filter.py

Filters the table to the organisms and samples the study is about, and reports what each
filter cost. Runs in the QIIME 2 environment. Standard library Python 3.8 or later.

```bash
python workflow/08_filter.py -b ~/research/baxter2016/q2/collapsed/table_by_sample.qza -r ~/research/baxter2016/q2/dada2/rep_seqs.qza -t ~/research/baxter2016/q2/taxonomy_final/taxonomy_gg2_2024.09_v4.qza -m ~/research/baxter2016/q2/collapsed/sample_metadata.tsv --drop-where "in_study_metadata='no'" --min-samples-fraction 0.05 -o ~/research/baxter2016/q2/filtered
```

Three filters in the order that makes each one interpretable: taxonomy first (keep the
target domain, drop mitochondria and chloroplast), then samples (mocks and controls out
before any community statistic is computed, because a mock left in shifts every
between-sample distance), then prevalence.

**The trap this stage was written around.** `qiime feature-table filter-features` has
`--p-filter-empty-samples` on by default, so filtering *features* can silently remove
*samples*. Every step reports samples before and after, and a sample lost to a feature
filter is called out with the reason.

**The prevalence threshold is a choice about what counts as evidence**, not a technical
step, so it lands in `filter_summary.tsv` next to the result along with every other
setting. The representative sequences are filtered against the final table, because a
tree built from sequences the table no longer holds is a silent mismatch that surfaces
much later as a diversity metric that cannot be computed.

| Option | Meaning |
|---|---|
| `--include`, `--exclude` | Taxa to keep and drop, default `Bacteria` and `mitochondria,chloroplast` |
| `--drop-where CLAUSE` | SQLite WHERE selecting samples to **drop**, e.g. `"in_study_metadata='no'"` |
| `--min-samples-fraction F` | Drop features seen in fewer than this fraction of samples, default 0, off |
| `--min-sample-reads N` | Drop samples below this many reads, default 0, off |

### Outputs

| File | Contents |
|---|---|
| `table_filtered.qza`, `rep_seqs_filtered.qza` | The filtered table and the sequences that match it |
| `filter_summary.tsv` | Every filter in order, what survived, and what it cost, with the settings at the top |
| `filter_per_sample.tsv` | Per sample: reads before, reads after, and whether it was dropped |
| `filter_log.txt` | Command, versions, each step and every sample or feature loss |

### Guards

| Guard | Why |
|---|---|
| Samples lost to a *feature* filter are named | `--p-filter-empty-samples` would otherwise change the sample count without a word |
| An emptied table is fatal, with the likely cause | An `--include` label that does not appear in the taxonomy removes everything |
| A `--drop-where` matching nothing warns | Usually a wrong column or value rather than a clean table |
| Sequences are filtered against the final table | Keeps the tree and the table in step |

---

## 09_tree.py

Aligns, masks, builds and roots the phylogeny, then proves its tips cover the table. Runs
in the QIIME 2 environment. Standard library Python 3.8 or later.

```bash
python workflow/09_tree.py -r ~/research/baxter2016/q2/filtered/rep_seqs_filtered.qza -b ~/research/baxter2016/q2/filtered/table_filtered.qza -o ~/research/baxter2016/q2/tree --threads 4
```

**The tips must match the table.** A phylogenetic diversity metric needs every feature in
the table to be a tip in the tree. Filter the table after building the tree and UniFrac
either fails much later with an opaque message or quietly computes on a subset, so this
stage compares the two sets and names the features that are missing.

**The masking step can remove most of the alignment.** On badly aligned input the masked
alignment collapses to a fraction of its length and the tree is built on almost nothing
with no error, so the length before and after masking is always reported and a floor is
applied.

| Option | Meaning |
|---|---|
| `-b`, `--table` | Feature table, to check the tips cover it. Omitting it is warned about loudly |
| `--threads N` | MAFFT and FastTree threads, default 1, always passed explicitly |
| `--mask-max-gap-frequency`, `--mask-min-conservation` | QIIME 2 defaults 1.0 and 0.4 |
| `--min-masked-fraction F` | Fail if masking leaves less than this much, default 0.25, our choice |

### Outputs

| File | Contents |
|---|---|
| `rooted_tree.qza`, `unrooted_tree.qza` | The trees, rooted at the midpoint and not |
| `alignment.qza`, `masked_alignment.qza` | The alignment before and after masking |
| `tree_tips.tsv` | Every id, and whether it is in the tree and in the table |
| `tree_log.txt` | Command, versions, alignment lengths and the tip comparison |

### Guards

| Guard | Why |
|---|---|
| Every table feature must be a tip | The reason the stage exists; a missing feature breaks UniFrac much later |
| A collapsed masked alignment is fatal | A tree built on it would be noise, and nothing raises an error |
| A sequence lost between alignment and masking is fatal | The two must hold the same records |
| Internal node labels are not counted as tips | FastTree writes support values where a label goes; counting them inflated 654 tips to 973 |

---

## 10_diversity.py

Alpha and beta diversity, with the rarefaction depth chosen from evidence. Runs in the
QIIME 2 environment. Standard library Python 3.8 or later.

```bash
python workflow/10_diversity.py -b ~/research/baxter2016/q2/filtered/table_filtered.qza -p ~/research/baxter2016/q2/tree/rooted_tree.qza -m ~/research/baxter2016/q2/collapsed/sample_metadata.tsv --coverage-table ~/research/baxter2016/q2/collapsed/table_by_sample.qza --depth 10000 --group-column dx -o ~/research/baxter2016/q2/diversity
```

**There is no published rule for a rarefaction depth.** Schloss 2024 (mSphere,
PMC10900887) found rarefaction the only approach that controls uneven sequencing effort
across common alpha and beta metrics, over datasets spanning 100-fold variation, and says
of the threshold: "My personal process for selecting a rarefaction threshold involves
looking for a natural break in the distribution of the number of sequences." So this
stage computes that break and tabulates what each candidate depth costs in samples.

**It refuses to invent a depth.** No `--depth` and no `--auto-depth` is an error.
`--auto-depth` uses the break, and if no break stands out it still refuses, because the
absence of a break is the finding: the choice is arbitrary there and the sensitivity of
any conclusion to it should be reported.

**`core-metrics-phylogenetic` subsamples once.** Schloss separates *rarefying*, a single
subsample, from *rarefaction*, repeating it 100 to 1,000 times and averaging, and argues
the conflation "was lost on many subsequent researchers". The single subsample is what
QIIME 2 gives you, so `alpha-rarefaction` runs alongside with `--p-iterations`, and both
the log and `diversity_settings.tsv` record which number came from which.

**Good's coverage is degenerate on ASV data, and the stage says so.** It is built on
singletons, and DADA2 does not emit singleton ASVs by design (benjjneb/dada2 issue 1491:
singletons are "too difficult to differentiate from errors"). So on a DADA2 table it
returns 1.0 for every sample, which looks like excellent depth and is a fact about the
denoiser. When the table holds no singletons the stage reports that rather than printing
a perfect score. Good's coverage belongs to OTU pipelines that keep singletons. A
prevalence filter has the same effect, which is why `--coverage-table` wants the
pre-filter table.

| Option | Meaning |
|---|---|
| `--depth N` | Rarefaction depth. Without it, nothing is guessed |
| `--auto-depth` | Use the natural break, and log which depth that gave |
| `--coverage-table FILE` | Table from **before** prevalence filtering, for Good's coverage |
| `--group-column NAME` | Metadata column to test groups on. Repeat for more than one |
| `--iterations N` | Rarefaction curve iterations, default 100 |
| `--max-sample-loss F` | Refuse a depth dropping more than this fraction, default 0.1, our choice |

### Outputs

| File | Contents |
|---|---|
| `core_metrics/` | The QIIME 2 alpha vectors, distance matrices and PCoA plots |
| `alpha_rarefaction.qzv` | The repeated-subsampling curves |
| `alpha_*.qzv`, `beta_*.qzv` | Group tests per metric and column, PERMANOVA with 999 permutations |
| `depth_candidates.tsv` | Each candidate depth, samples kept and lost, and reads used |
| `sample_depths.tsv` | Per sample: reads, whether it survives the depth, and Good's coverage |
| `diversity_settings.tsv` | The depth and its source, the iteration count, and the coverage status |
| `diversity_log.txt` | Command, versions, the break, the caveats and every group test |

### Guards

| Guard | Why |
|---|---|
| No depth is invented | There is no standard depth, and a default would be a fabricated one |
| A depth dropping more than `--max-sample-loss` is fatal | Discarding a tenth of a study should be deliberate and recorded |
| A depth above every sample is fatal | Nothing would be left |
| A table with no singletons is called degenerate | Good's coverage of 1.0 would otherwise read as evidence of depth |
| The single-subsample caveat is always logged | So a `core-metrics` number is never reported as a rarefied one |

---

## amplicon_regions.py

Not a stage. The region logic shared by `04_mock.py` and `06_resolution.py`: IUPAC-aware
primer matching, extracting what lies between a primer pair, and grouping sequences that
are identical over that region. One tested implementation rather than a copy in each
stage, so a fix to either reaches both.

Everything in it is IUPAC-aware, because primers carry ambiguity codes and treating them
as literal bases finds nothing while raising no error. `find_primer` prefers an exact
match anywhere over an earlier near match, since a near match exists to rescue references
whose primer site carries a real mismatch, not to license cutting at the first thing that
looks close.
