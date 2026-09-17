#!/usr/bin/env python3
"""Download raw FASTQ for one SRA/ENA accession and build a QIIME 2 manifest.

Give it one accession: a BioProject (PRJNA290926), a study (SRP062005), a
BioSample (SAMN03939374) or a run (SRR2144132). The script

  1. resolves the runs through the ENA Portal API,
  2. picks the read files for each run,
  3. downloads them over HTTPS with resume, checking byte size and MD5,
  4. pulls sample attributes from NCBI BioSample,
  5. writes a QIIME 2 manifest with one row per run,
  6. logs the command, versions, queries and checksums.

Standard library only. Needs Python 3.8 or later.

Example:
    python workflow/00_fetch_sra.py PRJNA290926 -o ~/research/baxter2016/raw
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import hashlib
import io
import logging
import os
import platform
import posixpath
import re
import shlex
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone

__version__ = "0.1.0"

ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ENA_FIELDS = [
    "run_accession", "sample_accession", "sample_title", "study_accession",
    "secondary_study_accession", "experiment_accession", "instrument_platform",
    "instrument_model", "library_layout", "library_strategy", "library_source",
    "read_count", "base_count", "fastq_ftp", "fastq_md5", "fastq_bytes",
    "first_public",
]
# One accession per query. ENA answers a comma-separated list, or an unknown
# ID, with HTTP 200 and an empty table, so both have to be caught here.
ACCESSION_RE = re.compile(r"^(PRJ(EB|NA|DB)\d+|[SED]R[APSXR]\d+|SAM(N|EA|D)\d+)$")
MD5_RE = re.compile(r"^[0-9a-f]{32}$")
NCBI_BATCH = 100
NCBI_MIN_INTERVAL = 0.34  # NCBI allows 3 requests per second without an API key
CHUNK = 1 << 20
BACKOFF = 2  # retry waits are BACKOFF ** attempt seconds, capped at 60
SRA_TIMEOUT = 3600  # seconds before fasterq-dump is killed
USER_AGENT = f"AmpliPub-fetch_sra/{__version__}"

log = logging.getLogger("fetch_sra")


class FetchError(RuntimeError):
    """A condition that must stop the run with a clear message."""


@dataclass
class ReadFile:
    run: str
    role: str  # "forward", "reverse" or "single"
    url: str
    md5: str
    size: int

    @property
    def name(self) -> str:
        return self.url.rsplit("/", 1)[-1]


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------

def _sleep(attempt: int) -> None:
    time.sleep(min(BACKOFF ** attempt, 60))


def http_get(url: str, retries: int = 3, timeout: int = 180) -> bytes:
    """GET a URL, retrying network errors. 4xx responses fail at once."""
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if 400 <= exc.code < 500:
                raise FetchError(f"HTTP {exc.code} for {url}") from exc
            err = exc
        except (urllib.error.URLError, OSError) as exc:
            err = exc
        if attempt == retries:
            raise FetchError(f"GET failed after {retries} attempts: {url}\n  {err}")
        log.warning("GET failed (%s), retrying: %s", err, url)
        _sleep(attempt)
    raise AssertionError("unreachable")


def _open(url: str, start: int):
    """Open a download stream from byte `start`. Returns (status, response)."""
    headers = {"User-Agent": USER_AGENT}
    if start:
        headers["Range"] = f"bytes={start}-"
    req = urllib.request.Request(url, headers=headers)
    resp = urllib.request.urlopen(req, timeout=120)
    return resp.status, resp


# --------------------------------------------------------------------------
# Stage 1: resolve runs
# --------------------------------------------------------------------------

def validate_accession(accession: str) -> str:
    acc = accession.strip()
    if "," in acc or " " in acc:
        raise FetchError(
            "Give one accession. ENA returns an empty table for a list, so use "
            "--runs to download a subset."
        )
    if not ACCESSION_RE.match(acc):
        raise FetchError(
            f"'{acc}' is not a BioProject, study, BioSample or run accession "
            "(expected e.g. PRJNA290926, SRP062005, SAMN03939374, SRR2144132)."
        )
    return acc


def ena_query_url(accession: str) -> str:
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
    }
    return ENA_FILEREPORT + "?" + urllib.parse.urlencode(params, safe=",")


def parse_runs(text: str, accession: str) -> list[dict]:
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise FetchError(
            f"ENA returned no runs for {accession}. ENA answers unknown or "
            "private accessions with an empty table, not an error, so check the "
            "ID and whether the data are public yet."
        )
    missing = [f for f in ENA_FIELDS if f not in rows[0]]
    if missing:
        raise FetchError(f"ENA response lacks expected columns: {', '.join(missing)}")
    dups = [r for r, n in Counter(r["run_accession"] for r in rows).items() if n > 1]
    if dups:
        raise FetchError(f"ENA listed runs more than once: {', '.join(dups[:10])}")
    return rows


def read_runs_file(path: str) -> list[str]:
    runs = []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if line:
                runs.append(line.split()[0])
    if not runs:
        raise FetchError(f"--runs file {path} lists no runs")
    return runs


def select_runs(rows: list[dict], wanted: list[str]) -> list[dict]:
    by_run = {r["run_accession"]: r for r in rows}
    unknown = [w for w in wanted if w not in by_run]
    if unknown:
        raise FetchError(
            f"{len(unknown)} run(s) in --runs are not part of this accession: "
            + ", ".join(unknown[:10])
        )
    return [by_run[w] for w in dict.fromkeys(wanted)]


def report_instruments(rows: list[dict]) -> None:
    """Log instrument labels with read lengths. Never filters on them."""
    stats: dict[str, list[int]] = {}
    for r in rows:
        s = stats.setdefault(r["instrument_model"] or "(blank)", [0, 0, 0])
        s[0] += 1
        s[1] += int(r["read_count"] or 0)
        s[2] += int(r["base_count"] or 0)
    log.info("Instrument labels as submitted (runs, mean bases per spot):")
    for model, (n, reads, bases) in sorted(stats.items()):
        mean = f"{bases / reads:.0f}" if reads else "NA"
        log.info("  %-24s %6d runs  %s", model, n, mean)
    if len(stats) > 1:
        log.warning(
            "Runs carry %d different instrument labels. SRA labels are sometimes "
            "wrong (Baxter 2016 has 467 Illumina runs labelled '454 GS'). Nothing "
            "is filtered on this field; compare the read lengths above.",
            len(stats),
        )


# --------------------------------------------------------------------------
# Stage 2: plan files
# --------------------------------------------------------------------------

def _split(value: str) -> list[str]:
    return [v for v in (value or "").split(";") if v]


def plan_files(row: dict) -> list[ReadFile]:
    """Choose the read files for one run. Empty list if ENA has none."""
    run = row["run_accession"]
    urls, md5s, sizes = (_split(row[k]) for k in ("fastq_ftp", "fastq_md5", "fastq_bytes"))
    if not urls:
        return []
    if not len(urls) == len(md5s) == len(sizes):
        raise FetchError(
            f"{run}: ENA lists {len(urls)} files but {len(md5s)} MD5s and {len(sizes)} sizes"
        )
    files: dict[str, ReadFile] = {}
    for url, md5, size in zip(urls, md5s, sizes):
        name = url.rsplit("/", 1)[-1]
        role = {f"{run}_1.fastq.gz": "forward", f"{run}_2.fastq.gz": "reverse",
                f"{run}.fastq.gz": "single"}.get(name)
        if role is None:
            raise FetchError(f"{run}: unexpected file name {name}")
        if not MD5_RE.match(md5):
            raise FetchError(f"{run}: bad MD5 '{md5}' for {name}")
        if not url.startswith("http"):
            url = "https://" + url
        files[role] = ReadFile(run, role, url, md5, int(size))
    if "forward" in files or "reverse" in files:
        if not ("forward" in files and "reverse" in files):
            present = "forward" if "forward" in files else "reverse"
            raise FetchError(
                f"{run}: only the {present} file is listed. Refusing to treat a "
                "broken pair as single-end."
            )
        if "single" in files:
            log.info("%s: skipping unpaired file %s (%d bytes); using _1 and _2",
                     run, files["single"].name, files["single"].size)
        return [files["forward"], files["reverse"]]
    return [files["single"]]


def check_layout(plan: dict[str, list[ReadFile]], rows: list[dict]) -> str:
    layouts = {run: ("paired" if len(f) == 2 else "single") for run, f in plan.items() if f}
    kinds = Counter(layouts.values())
    if len(kinds) > 1:
        raise FetchError(
            f"Mixed layouts: {kinds['paired']} paired and {kinds['single']} single-end "
            "runs. Use --runs to download one layout at a time."
        )
    layout = next(iter(kinds))
    for r in rows:
        declared = (r["library_layout"] or "").lower()
        actual = layouts.get(r["run_accession"])
        if actual and declared and declared != actual:
            log.warning("%s: ENA says %s but the files are %s; using the files",
                        r["run_accession"], declared, actual)
    return layout


# --------------------------------------------------------------------------
# Stage 3: download
# --------------------------------------------------------------------------

def md5sum(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def download(rf: ReadFile, dest_dir: str, retries: int = 5) -> str:
    """Fetch one file with resume. Returns 'present' or 'downloaded'.

    Nothing reaches its final name until size and MD5 both match ENA.
    """
    final = os.path.join(dest_dir, rf.name)
    if os.path.exists(final):
        if os.path.getsize(final) == rf.size and md5sum(final) == rf.md5:
            return "present"
        log.warning("%s exists but fails size or MD5; downloading again", rf.name)
        os.remove(final)
    part = final + ".part"
    for attempt in range(1, retries + 1):
        start = os.path.getsize(part) if os.path.exists(part) else 0
        if start > rf.size:
            os.remove(part)
            start = 0
        if start < rf.size:
            if start:
                log.info("%s: resuming at byte %d of %d", rf.name, start, rf.size)
            try:
                status, resp = _open(rf.url, start)
                with resp:
                    if start and status != 206:  # server ignored the range
                        start = 0
                    with open(part, "ab" if start else "wb") as out:
                        shutil.copyfileobj(resp, out, CHUNK)
            except urllib.error.HTTPError as exc:
                if exc.code == 416:  # range not satisfiable, start over
                    os.remove(part)
                elif 400 <= exc.code < 500:
                    raise FetchError(f"{rf.name}: HTTP {exc.code} for {rf.url}") from exc
                log.warning("%s: HTTP %s (attempt %d of %d)", rf.name, exc.code, attempt, retries)
                _sleep(attempt)
                continue
            except (urllib.error.URLError, OSError) as exc:
                log.warning("%s: %s (attempt %d of %d)", rf.name, exc, attempt, retries)
                _sleep(attempt)
                continue
        size = os.path.getsize(part)
        if size != rf.size:
            log.warning("%s: have %d of %d bytes (attempt %d of %d)",
                        rf.name, size, rf.size, attempt, retries)
            if size > rf.size:
                os.remove(part)
            _sleep(attempt)
            continue
        if md5sum(part) != rf.md5:
            log.warning("%s: MD5 mismatch, discarding (attempt %d of %d)", rf.name, attempt, retries)
            os.remove(part)
            continue
        os.replace(part, final)
        return "downloaded"
    if os.path.exists(part) and os.path.getsize(part) >= rf.size:
        os.remove(part)  # complete but wrong; a shorter .part is kept for resume
    raise FetchError(f"{rf.name}: failed size/MD5 verification after {retries} attempts")


def download_all(files: list[ReadFile], dest_dir: str, threads: int, retries: int,
                 collect_failures: bool = False) -> list[ReadFile]:
    """Download every file. Returns the files that failed verification.

    With collect_failures False one failure is fatal. With it True the caller
    decides what to do, which is how --sra-fallback reaches files ENA lists in
    its filereport but will not actually serve.
    """
    failures: list[str] = []
    failed: list[ReadFile] = []
    done = 0
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = {pool.submit(download, rf, dest_dir, retries): rf for rf in files}
        for fut in as_completed(futures):
            rf = futures[fut]
            done += 1
            try:
                status = fut.result()
                log.info("[%d/%d] %s %s", done, len(files), status, rf.name)
            except FetchError as exc:
                failures.append(str(exc))
                failed.append(rf)
                log.error("[%d/%d] FAILED %s", done, len(files), rf.name)
    if failures and not collect_failures:
        raise FetchError(f"{len(failures)} file(s) failed:\n  " + "\n  ".join(failures)
                         + "\nRerun the same command to resume. If ENA lists a file "
                         "but will not serve it, rerun with --sra-fallback.")
    return failed


def fetch_with_sra_tools(row: dict, dest_dir: str, threads: int) -> list[ReadFile]:
    """Fetch a whole run from NCBI SRA.

    Used for runs ENA does not mirror and for runs whose ENA files fail
    verification. The whole run is refetched, not the single broken mate, so
    both mates come from one source and stay in the same read order. MD5s
    here are computed locally, so the check is a read count against ENA
    rather than a checksum against the archive."""
    run = row["run_accession"]
    exe = shutil.which("fasterq-dump")
    if not exe:
        raise FetchError(f"{run}: ENA has no FASTQ and fasterq-dump is not on PATH")
    tmp = os.path.join(dest_dir, f".tmp_{run}")
    os.makedirs(tmp, exist_ok=True)
    for stale in glob.glob(os.path.join(dest_dir, f"{run}*.fastq.gz.part")):
        os.remove(stale)  # a half-finished ENA download of the same run
    cmd = [exe, "--split-files", "--threads", str(threads), "--temp", tmp, "--outdir", tmp, run]
    log.info("%s: %s", run, shlex.join(cmd))
    try:
        # stdin closed and a time limit, so a stalled download cannot hang the run.
        subprocess.run(cmd, check=True, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, timeout=SRA_TIMEOUT)
    except subprocess.TimeoutExpired as exc:
        raise FetchError(f"{run}: fasterq-dump timed out after {SRA_TIMEOUT} s and was killed") from exc
    except subprocess.CalledProcessError as exc:
        raise FetchError(f"{run}: fasterq-dump failed:\n{exc.stdout.decode(errors='replace')}") from exc
    produced = sorted(f for f in os.listdir(tmp) if f.endswith(".fastq"))
    roles = {f"{run}_1.fastq": "forward", f"{run}_2.fastq": "reverse", f"{run}.fastq": "single"}
    unknown = [f for f in produced if f not in roles]
    if unknown:
        raise FetchError(f"{run}: fasterq-dump wrote unexpected files: {unknown}")
    if f"{run}_1.fastq" in produced:
        produced = [f"{run}_1.fastq", f"{run}_2.fastq"]
        if not os.path.exists(os.path.join(tmp, produced[1])):
            raise FetchError(f"{run}: fasterq-dump wrote _1 without _2")
    out = []
    for name in produced:
        src = os.path.join(tmp, name)
        with open(src, "rb") as fh:
            n_reads = sum(1 for _ in fh) // 4
        expected = int(row["read_count"] or 0)
        if expected and n_reads != expected:
            raise FetchError(f"{run}: {name} has {n_reads} reads, ENA says {expected}")
        dst = os.path.join(dest_dir, name + ".gz")
        with open(src, "rb") as fin, gzip.open(dst + ".part", "wb") as fout:
            shutil.copyfileobj(fin, fout, CHUNK)
        os.replace(dst + ".part", dst)
        out.append(ReadFile(run, roles[name], "file://" + dst, md5sum(dst), os.path.getsize(dst)))
    shutil.rmtree(tmp)
    return out


# --------------------------------------------------------------------------
# Stage 4: sample metadata
# --------------------------------------------------------------------------

def parse_biosample_xml(data: bytes) -> dict[str, dict[str, str]]:
    samples = {}
    for bs in ET.fromstring(data).iter("BioSample"):
        attrs = {"biosample_title": (bs.findtext("Description/Title") or "").strip()}
        for a in bs.iter("Attribute"):
            if a.get("attribute_name"):
                attrs[a.get("attribute_name")] = (a.text or "").strip()
        samples[bs.get("accession")] = attrs
    return samples


def fetch_biosamples(accessions: list[str]) -> dict[str, dict[str, str]]:
    accs = sorted({a for a in accessions if a})
    out: dict[str, dict[str, str]] = {}
    last = 0.0
    for i in range(0, len(accs), NCBI_BATCH):
        wait = NCBI_MIN_INTERVAL - (time.monotonic() - last)
        if wait > 0:
            time.sleep(wait)
        last = time.monotonic()
        params = {"db": "biosample", "id": ",".join(accs[i:i + NCBI_BATCH]), "retmode": "xml"}
        out.update(parse_biosample_xml(http_get(NCBI_EFETCH + "?" + urllib.parse.urlencode(params, safe=","))))
    return out


def _clean(value: str) -> str:
    return re.sub(r"[\t\r\n]+", " ", value or "").strip()


def write_sample_metadata(path: str, rows: list[dict], biosamples: dict) -> None:
    samples: dict[str, dict] = {}
    for r in rows:
        s = samples.setdefault(r["sample_accession"], {"sample_title": r["sample_title"], "runs": []})
        s["runs"].append(r["run_accession"])
    missing = [s for s in samples if s not in biosamples]
    if missing:
        log.warning("%d of %d samples had no BioSample record (e.g. %s); their rows "
                    "carry ENA fields only", len(missing), len(samples), missing[0])
    titles = Counter(s["sample_title"] for s in samples.values() if s["sample_title"])
    shared = [t for t, n in titles.items() if n > 1]
    if shared:
        log.warning("%d sample titles are shared by different BioSamples (e.g. '%s'); "
                    "group by sample accession, not title", len(shared), shared[0])
    attr_names: list[str] = []
    for acc in samples:
        for k in biosamples.get(acc, {}):
            if k not in attr_names:
                attr_names.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "sample_title", "runs"] + attr_names)
        for acc, s in samples.items():
            bs = biosamples.get(acc, {})
            w.writerow([acc, _clean(s["sample_title"]), ";".join(s["runs"])]
                       + [_clean(bs.get(k, "")) for k in attr_names])


# --------------------------------------------------------------------------
# Stage 5 and 6: manifest, map, checksums, log
# --------------------------------------------------------------------------

MANIFEST_HEADER = {
    # Column names from q2_types.per_sample_sequences._formats, QIIME 2 2025.7
    "paired": ["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"],
    "single": ["sample-id", "absolute-filepath"],
}


def write_manifest(path: str, plan: dict[str, list[ReadFile]], layout: str, fastq_dir: str) -> None:
    """Paths are written as $PWD/<fastq dir>/<file>, relative to the manifest's directory.

    An absolute path ties the download to the machine it was made on; moving it to another
    machine or into a container broke the import (2026-09-17). QIIME 2 2025.7 rejects plain
    relative paths ("must be absolute") but expands $PWD, so the manifest works from any
    location as long as the import runs with the manifest's directory as the working
    directory, which the Snakefile's import rule does.
    """
    rel = os.path.relpath(os.path.abspath(fastq_dir), os.path.dirname(os.path.abspath(path)))
    if rel == os.pardir or rel.startswith(os.pardir + os.sep) or os.path.isabs(rel):
        raise FetchError(f"{fastq_dir} is not inside the manifest's directory, "
                         "so the manifest could not be moved with it")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(MANIFEST_HEADER[layout])
        for run, files in plan.items():
            by_role = {f.role: "$PWD/" + posixpath.join(*rel.split(os.sep), f.name)
                       for f in files}
            if layout == "paired":
                w.writerow([run, by_role["forward"], by_role["reverse"]])
            else:
                w.writerow([run, by_role["single"]])


def write_run_map(path: str, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample-id", "sample_accession", "sample_title", "instrument_model", "read_count"])
        for r in rows:
            w.writerow([r["run_accession"], r["sample_accession"], _clean(r["sample_title"]),
                        r["instrument_model"], r["read_count"]])


def write_run_sources(path: str, rows: list[dict], from_sra: set[str]) -> None:
    """Record where each run's FASTQ actually came from, for the methods section."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "source", "verified_by"])
        for r in rows:
            run = r["run_accession"]
            if run in from_sra:
                w.writerow([run, "sra", "read_count"])
            else:
                w.writerow([run, "ena", "size+md5"])


def write_rows(path: str, rows: list[dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ENA_FIELDS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "fetch_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Download raw FASTQ for one SRA/ENA accession and build a QIIME 2 manifest.")
    p.add_argument("accession", help="BioProject, study, BioSample or run (e.g. PRJNA290926)")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("--runs", help="file of run accessions to fetch, one per line (subset)")
    p.add_argument("--expect-runs", type=int,
                   help="fail unless the accession resolves to exactly this many runs")
    p.add_argument("--threads", type=int, default=4, help="parallel downloads (default 4)")
    p.add_argument("--retries", type=int, default=5, help="attempts per file (default 5)")
    p.add_argument("--dry-run", action="store_true", help="resolve and plan, download nothing")
    p.add_argument("--sra-fallback", action="store_true",
                   help="use fasterq-dump for runs ENA does not mirror, and for runs "
                        "whose ENA files fail size or MD5 verification")
    p.add_argument("--no-metadata", action="store_true", help="skip NCBI BioSample attributes")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    accession = validate_accession(args.accession)
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    fastq_dir = os.path.join(outdir, "fastq")
    os.makedirs(fastq_dir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("fetch_sra %s | Python %s | %s", __version__, platform.python_version(), platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    # Stage 1
    url = ena_query_url(accession)
    log.info("ENA query: %s", url)
    text = http_get(url).decode("utf-8")
    with open(os.path.join(outdir, "ena_filereport.tsv"), "w", newline="") as fh:
        fh.write(text)
    rows = parse_runs(text, accession)
    log.info("%s resolves to %d runs across %d samples", accession, len(rows),
             len({r["sample_accession"] for r in rows}))
    if args.expect_runs is not None and len(rows) != args.expect_runs:
        raise FetchError(f"expected {args.expect_runs} runs, ENA lists {len(rows)}")
    if args.runs:
        rows = select_runs(rows, read_runs_file(args.runs))
        log.info("--runs keeps %d runs across %d samples", len(rows),
                 len({r["sample_accession"] for r in rows}))
    report_instruments(rows)
    write_rows(os.path.join(outdir, "runs.tsv"), rows)

    # Stage 2
    plan = {r["run_accession"]: plan_files(r) for r in rows}
    no_fastq = [run for run, f in plan.items() if not f]
    if no_fastq and not args.sra_fallback:
        raise FetchError(f"{len(no_fastq)} run(s) have no FASTQ on ENA (e.g. {no_fastq[0]}). "
                         "Rerun with --sra-fallback to fetch them with fasterq-dump.")
    layout = check_layout(plan, rows) if any(plan.values()) else "paired"
    files = [f for group in plan.values() for f in group]
    total = sum(f.size for f in files)
    log.info("plan: %d %s-end runs, %d files, %.2f GB from ENA; %d run(s) via sra-tools",
             len(rows) - len(no_fastq), layout, len(files), total / 1e9, len(no_fastq))
    if args.dry_run:
        log.info("dry run: nothing downloaded. Plan written to runs.tsv")
        return

    # Stage 3
    failed = download_all(files, fastq_dir, args.threads, args.retries,
                          collect_failures=args.sra_fallback)
    by_row = {r["run_accession"]: r for r in rows}
    broken = list(dict.fromkeys(f.run for f in failed))
    if broken:
        log.warning("%d run(s) have files ENA lists but will not serve correctly: %s. "
                    "Refetching each whole run from SRA so both mates come from one "
                    "source and keep the same read order.", len(broken), ", ".join(broken))
    from_sra = []
    for run_acc in no_fastq + broken:
        plan[run_acc] = fetch_with_sra_tools(by_row[run_acc], fastq_dir, args.threads)
        from_sra.append(run_acc)
    files = [f for group in plan.values() for f in group]
    layout = check_layout(plan, rows)

    # Stage 4
    if args.no_metadata:
        biosamples = {}
    else:
        log.info("fetching BioSample attributes from NCBI")
        biosamples = fetch_biosamples([r["sample_accession"] for r in rows])
    write_sample_metadata(os.path.join(outdir, "sample_metadata.tsv"), rows, biosamples)

    # Stage 5 and 6
    manifest = os.path.join(outdir, "manifest.tsv")
    write_manifest(manifest, plan, layout, fastq_dir)
    write_run_map(os.path.join(outdir, "run_to_sample.tsv"), rows)
    write_run_sources(os.path.join(outdir, "run_sources.tsv"), rows, set(from_sra))
    with open(os.path.join(outdir, "checksums.md5"), "w", newline="") as fh:
        for f in sorted(files, key=lambda f: f.name):
            fh.write(f"{f.md5}  fastq/{f.name}\n")
    if from_sra:
        log.warning("%d run(s) came from SRA, not ENA: %s. Their MD5s in "
                    "checksums.md5 were computed here, so they verify the local copy "
                    "only. The archive check for them was the read count. See "
                    "run_sources.tsv.", len(from_sra), ", ".join(sorted(from_sra)))
    n_rows = sum(1 for _ in open(manifest)) - 1
    if n_rows != len(rows):
        raise FetchError(f"manifest has {n_rows} rows for {len(rows)} runs")
    log.info("done: %d runs, %d files, %.2f GB. Verify with: cd %s && md5sum -c checksums.md5",
             len(rows), len(files), sum(f.size for f in files) / 1e9, shlex.quote(outdir))
    log.info("import (the manifest uses $PWD, so run it from %s): cd %s && qiime tools import "
             "--type 'SampleData[%s]' --input-format %s --input-path manifest.tsv "
             "--output-path demux.qza",
             outdir, shlex.quote(outdir),
             "PairedEndSequencesWithQuality" if layout == "paired" else "SequencesWithQuality",
             "PairedEndFastqManifestPhred33V2" if layout == "paired" else "SingleEndFastqManifestPhred33V2")


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except FetchError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
