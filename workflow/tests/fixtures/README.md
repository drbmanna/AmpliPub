# Test fixtures

Real API responses saved on 2026-09-11. Not edited, only trimmed to a few rows.

| File | Source | Why it is here |
|---|---|---|
| `ena_header_only.tsv` | ENA filereport for `PRJNA000000000` | ENA answers an unknown accession with HTTP 200 and a header-only table |
| `ena_baxter_subset.tsv` | ENA filereport for `PRJNA290926`, rows SRR2143519, SRR2143538, SRR2144132, SRR2143955, SRR2143956 | One run mislabelled `454 GS`, one labelled `Illumina MiSeq`, one mock community, and both runs of a resequenced sample |
| `ena_three_files.tsv` | ENA filereport for `DRR045356` | A paired run that also lists an unpaired `DRR045356.fastq.gz` |
| `biosample_two.xml` | NCBI efetch, `db=biosample`, `SAMN03939374,SAMN03939742` | Submitter attributes for two Baxter 2016 samples |
| `SRR2144126_1_fastqc.zip` | FastQC 0.12.1 on `SRR2144126_1.fastq.gz`, run 2026-09-12 | Real report layout for `00_qc_raw.py`. Only `summary.txt` and `fastqc_data.txt` are kept, unedited; images and HTML are dropped for size |

The ENA field list matches `ENA_FIELDS` in `00_fetch_sra.py`.
