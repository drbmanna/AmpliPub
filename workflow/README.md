# Workflow scripts

Staged scripts that run upstream of the AmpliPub R package, in WSL or Linux. They are
not part of the package build.

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
