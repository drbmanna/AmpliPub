"""Offline tests for 00_fetch_sra.py. Each guard is shown to fire.

Fixtures are real ENA and NCBI responses saved on 2026-09-11 (see
fixtures/README.md). Downloads are served from memory, so no network is used.
"""

import csv
import hashlib
import importlib.util
import io
import logging
import os
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
FIX = HERE / "fixtures"
spec = importlib.util.spec_from_file_location("fetch_sra", HERE.parent / "00_fetch_sra.py")
fs = importlib.util.module_from_spec(spec)
sys.modules["fetch_sra"] = fs
spec.loader.exec_module(fs)
fs.BACKOFF = 0


def rows_from(name):
    return fs.parse_runs((FIX / name).read_text(), name)


def by_run(rows):
    return {r["run_accession"]: r for r in rows}


class FakeServer:
    """Stands in for fs._open. Serves bytes, honours Range unless told not to."""

    def __init__(self, blobs, truncate_first=False, always_truncate=False, ignore_range=False):
        self.blobs = blobs
        self.truncate_first = truncate_first
        self.always_truncate = always_truncate
        self.ignore_range = ignore_range
        self.calls = []

    def __call__(self, url, start):
        self.calls.append((url, start))
        eff = 0 if self.ignore_range else start
        body = self.blobs[url][eff:]
        if self.always_truncate or (self.truncate_first and len(self.calls) == 1):
            body = body[: len(body) // 2]
        return (206 if eff else 200), io.BytesIO(body)


def readfile(payload, name="SRRX_1.fastq.gz", md5=None):
    return fs.ReadFile("SRRX", "forward", "https://host/" + name,
                       md5 or hashlib.md5(payload).hexdigest(), len(payload))


# ---- accession and ENA response guards ----------------------------------

@pytest.mark.parametrize("acc", ["PRJNA290926", "SRP062005", "SAMN03939374", "SRR2144132",
                                 "ERR123", "DRR045356", "PRJEB1234", "SAMEA123"])
def test_valid_accessions_pass(acc):
    assert fs.validate_accession(acc) == acc


@pytest.mark.parametrize("acc", ["SRR1,SRR2", "SRR1 SRR2", "GSE12345", "prjna290926", ""])
def test_bad_accessions_fire(acc):
    with pytest.raises(fs.FetchError):
        fs.validate_accession(acc)


def test_header_only_response_is_fatal():
    with pytest.raises(fs.FetchError, match="no runs"):
        rows_from("ena_header_only.tsv")


def test_missing_columns_fatal():
    with pytest.raises(fs.FetchError, match="lacks expected columns"):
        fs.parse_runs("run_accession\tfastq_ftp\nSRR1\tx\n", "X")


def test_duplicate_runs_fatal():
    text = (FIX / "ena_baxter_subset.tsv").read_text()
    line = text.splitlines()[1]
    with pytest.raises(fs.FetchError, match="more than once"):
        fs.parse_runs(text + line + "\n", "X")


def test_unknown_run_in_subset_fatal():
    with pytest.raises(fs.FetchError, match="not part of this accession"):
        fs.select_runs(rows_from("ena_baxter_subset.tsv"), ["SRR2143519", "SRR9999999"])


def test_subset_keeps_order_and_drops_repeats():
    got = fs.select_runs(rows_from("ena_baxter_subset.tsv"), ["SRR2144132", "SRR2143519", "SRR2144132"])
    assert [r["run_accession"] for r in got] == ["SRR2144132", "SRR2143519"]


def test_runs_file_ignores_comments(tmp_path):
    p = tmp_path / "runs.txt"
    p.write_text("# dev subset\nSRR1  mock\n\nSRR2 # note\n")
    assert fs.read_runs_file(str(p)) == ["SRR1", "SRR2"]


# ---- file planning -------------------------------------------------------

def test_mislabelled_instrument_is_kept_and_warned(caplog):
    rows = rows_from("ena_baxter_subset.tsv")
    r = by_run(rows)["SRR2143519"]
    assert r["instrument_model"] == "454 GS"  # the real SRA label
    with caplog.at_level(logging.INFO, logger="fetch_sra"):
        fs.report_instruments(rows)
    assert "different instrument labels" in caplog.text
    assert [f.role for f in fs.plan_files(r)] == ["forward", "reverse"]


def test_three_file_run_keeps_only_mates():
    r = rows_from("ena_three_files.tsv")[0]
    files = fs.plan_files(r)
    assert [f.name for f in files] == ["DRR045356_1.fastq.gz", "DRR045356_2.fastq.gz"]
    assert all(f.url.startswith("https://ftp.sra.ebi.ac.uk/") for f in files)


def test_missing_mate_fatal():
    r = dict(by_run(rows_from("ena_baxter_subset.tsv"))["SRR2143519"])
    for k in ("fastq_ftp", "fastq_md5", "fastq_bytes"):
        r[k] = r[k].split(";")[0]
    with pytest.raises(fs.FetchError, match="broken pair"):
        fs.plan_files(r)


def test_file_md5_size_count_mismatch_fatal():
    r = dict(by_run(rows_from("ena_baxter_subset.tsv"))["SRR2143519"])
    r["fastq_md5"] = r["fastq_md5"].split(";")[0]
    with pytest.raises(fs.FetchError, match="2 files but 1 MD5s"):
        fs.plan_files(r)


def test_no_fastq_gives_empty_plan():
    r = dict(by_run(rows_from("ena_baxter_subset.tsv"))["SRR2143519"])
    r.update(fastq_ftp="", fastq_md5="", fastq_bytes="")
    assert fs.plan_files(r) == []


def test_mixed_layout_fatal():
    rows = rows_from("ena_baxter_subset.tsv")
    plan = {r["run_accession"]: fs.plan_files(r) for r in rows}
    first = next(iter(plan))
    plan[first] = plan[first][:1]
    with pytest.raises(fs.FetchError, match="Mixed layouts"):
        fs.check_layout(plan, rows)


# ---- download ------------------------------------------------------------

def test_download_then_rerun_skips(tmp_path, monkeypatch):
    payload = b"@r\nACGT\n+\nIIII\n" * 500
    rf = readfile(payload)
    server = FakeServer({rf.url: payload})
    monkeypatch.setattr(fs, "_open", server)
    assert fs.download(rf, str(tmp_path)) == "downloaded"
    assert fs.download(rf, str(tmp_path)) == "present"
    assert len(server.calls) == 1
    assert (tmp_path / rf.name).read_bytes() == payload


def test_md5_mismatch_fires_and_leaves_nothing(tmp_path, monkeypatch):
    payload = b"x" * 1000
    rf = readfile(payload, md5="0" * 32)
    monkeypatch.setattr(fs, "_open", FakeServer({rf.url: payload}))
    with pytest.raises(fs.FetchError, match="size/MD5"):
        fs.download(rf, str(tmp_path), retries=3)
    assert os.listdir(tmp_path) == []


def test_truncated_download_resumes_from_offset(tmp_path, monkeypatch):
    payload = bytes(range(256)) * 40
    rf = readfile(payload)
    server = FakeServer({rf.url: payload}, truncate_first=True)
    monkeypatch.setattr(fs, "_open", server)
    assert fs.download(rf, str(tmp_path)) == "downloaded"
    assert [s for _, s in server.calls] == [0, len(payload) // 2]
    assert (tmp_path / rf.name).read_bytes() == payload


def test_persistent_truncation_fires(tmp_path, monkeypatch):
    payload = b"y" * 1000
    rf = readfile(payload)
    monkeypatch.setattr(fs, "_open", FakeServer({rf.url: payload}, always_truncate=True))
    with pytest.raises(fs.FetchError, match="size/MD5"):
        fs.download(rf, str(tmp_path), retries=3)
    assert not (tmp_path / rf.name).exists()


def test_server_ignoring_range_restarts_cleanly(tmp_path, monkeypatch):
    payload = bytes(range(256)) * 10
    rf = readfile(payload)
    (tmp_path / (rf.name + ".part")).write_bytes(payload[:100])
    monkeypatch.setattr(fs, "_open", FakeServer({rf.url: payload}, ignore_range=True))
    assert fs.download(rf, str(tmp_path)) == "downloaded"
    assert (tmp_path / rf.name).read_bytes() == payload


def test_corrupt_existing_file_is_replaced(tmp_path, monkeypatch):
    payload = b"z" * 500
    rf = readfile(payload)
    (tmp_path / rf.name).write_bytes(b"w" * 500)  # right size, wrong content
    monkeypatch.setattr(fs, "_open", FakeServer({rf.url: payload}))
    assert fs.download(rf, str(tmp_path)) == "downloaded"
    assert (tmp_path / rf.name).read_bytes() == payload


# ---- metadata and outputs ------------------------------------------------

def test_biosample_attributes_parse():
    bs = fs.parse_biosample_xml((FIX / "biosample_two.xml").read_bytes())
    assert set(bs) == {"SAMN03939374", "SAMN03939742"}
    assert bs["SAMN03939374"]["diagnosis"] == "Adenoma"
    assert bs["SAMN03939374"]["biosample_title"] == "2009650"


def test_manifest_one_row_per_run_absolute_paths(tmp_path):
    rows = rows_from("ena_baxter_subset.tsv")
    plan = {r["run_accession"]: fs.plan_files(r) for r in rows}
    path = tmp_path / "manifest.tsv"
    fs.write_manifest(str(path), plan, "paired", str(tmp_path / "fastq"))
    table = list(csv.reader(path.open(), delimiter="\t"))
    assert table[0] == ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"]
    assert [t[0] for t in table[1:]] == list(plan)
    assert {"SRR2143955", "SRR2143956"} <= {t[0] for t in table[1:]}  # resequenced sample
    assert all(os.path.isabs(p) for t in table[1:] for p in t[1:])


def test_sample_metadata_groups_runs_and_warns_missing(tmp_path, caplog):
    rows = rows_from("ena_baxter_subset.tsv")
    bs = fs.parse_biosample_xml((FIX / "biosample_two.xml").read_bytes())
    path = tmp_path / "sample_metadata.tsv"
    with caplog.at_level(logging.WARNING, logger="fetch_sra"):
        fs.write_sample_metadata(str(path), rows, bs)
    assert "had no BioSample record" in caplog.text
    table = {t["sample-id"]: t for t in csv.DictReader(path.open(), delimiter="\t")}
    assert table["SAMN03939742"]["runs"].split(";") == ["SRR2143955", "SRR2143956"]
    assert table["SAMN03939374"]["diagnosis"] == "Adenoma"


# ---- end to end, offline -------------------------------------------------

def synthetic_project():
    """Two runs of one sample plus one run of another, payloads in memory."""
    blobs, lines = {}, ["\t".join(fs.ENA_FIELDS)]
    for run, sample in [("SRR1", "SAMN1"), ("SRR2", "SAMN1"), ("SRR3", "SAMN2")]:
        paths, md5s, sizes = [], [], []
        for mate in (1, 2):
            data = f"@{run}.{mate}\nACGT\n+\nIIII\n".encode() * 50
            path = f"ftp.sra.ebi.ac.uk/vol1/fastq/{run}_{mate}.fastq.gz"
            blobs["https://" + path] = data
            paths.append(path)
            md5s.append(hashlib.md5(data).hexdigest())
            sizes.append(str(len(data)))
        row = dict.fromkeys(fs.ENA_FIELDS, "")
        row.update(run_accession=run, sample_accession=sample, sample_title=sample.lower(),
                   instrument_model="Illumina MiSeq", library_layout="PAIRED",
                   read_count="50", base_count="400", fastq_ftp=";".join(paths),
                   fastq_md5=";".join(md5s), fastq_bytes=";".join(sizes))
        lines.append("\t".join(row[f] for f in fs.ENA_FIELDS))
    return "\n".join(lines) + "\n", blobs


def fake_http(ena_text):
    xml = (b'<BioSampleSet><BioSample accession="SAMN1"><Description><Title>s1</Title>'
           b'</Description><Attributes><Attribute attribute_name="dx">case</Attribute>'
           b'</Attributes></BioSample></BioSampleSet>')
    return lambda url, **kw: ena_text.encode() if "ebi.ac.uk" in url else xml


def test_end_to_end_offline(tmp_path, monkeypatch):
    ena_text, blobs = synthetic_project()
    monkeypatch.setattr(fs, "http_get", fake_http(ena_text))
    monkeypatch.setattr(fs, "_open", FakeServer(blobs))
    assert fs.main(["PRJNA1", "-o", str(tmp_path), "--expect-runs", "3"]) == 0
    manifest = list(csv.reader((tmp_path / "manifest.tsv").open(), delimiter="\t"))
    assert [m[0] for m in manifest[1:]] == ["SRR1", "SRR2", "SRR3"]
    sums = (tmp_path / "checksums.md5").read_text().splitlines()
    assert len(sums) == 6
    for line in sums:
        md5, rel = line.split("  ")
        assert hashlib.md5((tmp_path / rel).read_bytes()).hexdigest() == md5
    meta = {t["sample-id"]: t for t in csv.DictReader((tmp_path / "sample_metadata.tsv").open(), delimiter="\t")}
    assert meta["SAMN1"]["runs"] == "SRR1;SRR2" and meta["SAMN1"]["dx"] == "case"
    assert "command:" in (tmp_path / "fetch_log.txt").read_text()


def test_expect_runs_mismatch_fails(tmp_path, monkeypatch):
    ena_text, _ = synthetic_project()
    monkeypatch.setattr(fs, "http_get", fake_http(ena_text))
    assert fs.main(["PRJNA1", "-o", str(tmp_path), "--expect-runs", "4", "--dry-run"]) == 1


def test_bogus_accession_fails_through_main(tmp_path, monkeypatch):
    header_only = (FIX / "ena_header_only.tsv").read_text()
    monkeypatch.setattr(fs, "http_get", lambda url, **kw: header_only.encode())
    assert fs.main(["PRJNA000000000", "-o", str(tmp_path)]) == 1


def test_no_fastq_without_fallback_fails(tmp_path, monkeypatch):
    ena_text, _ = synthetic_project()
    lines = ena_text.splitlines()
    cols = lines[1].split("\t")
    for f in ("fastq_ftp", "fastq_md5", "fastq_bytes"):
        cols[fs.ENA_FIELDS.index(f)] = ""
    lines[1] = "\t".join(cols)
    monkeypatch.setattr(fs, "http_get", fake_http("\n".join(lines) + "\n"))
    assert fs.main(["PRJNA1", "-o", str(tmp_path), "--dry-run"]) == 1
