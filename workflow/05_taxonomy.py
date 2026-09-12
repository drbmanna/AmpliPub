#!/usr/bin/env python3
"""Classify ASVs against a reference taxonomy and report how deep the names go.

Give it the sequences from 03_dada2.py and a trained classifier. The script

  1. checks each classifier really is a TaxonomicClassifier before spending compute,
  2. runs qiime feature-classifier classify-sklearn for each one,
  3. reports how deep each one classifies, by ASV and by read,
  4. optionally, if given a second classifier, reports where their labels differ,
  5. optionally lists what each one called the mock community's exact-match ASVs.

**One reference is the default, and it should be.** Pick a reference, pin its version,
name it in the methods, and report the coverage table: the fraction of ASVs, and of
reads, that receive a name at each rank. That is the number a reader needs, and it is
what this stage produces when you pass a single classifier.

**A second classifier is a diagnostic, not a better default.** Comparing two references
tells you how much of your naming depends on the reference you chose, which is worth
knowing once. It is not a routine part of an analysis, and the comparison measures less
than it appears to: see `label_differences` below. Running two references as standard
also reintroduces the nomenclature problem that choosing one removes at the source.

**The ceiling is the fragment, not the classifier.** 253 bp of V4 does not carry
species-level information for many taxa. The reference work for this project found that
*S. aureus.1* and *S. epidermidis.1/2* have identical V4 sequences, so no method
separates them here. Read a low species-level rate as a property of the amplicon and
report it that way, rather than reaching for a cleverer model.

**Mock calls are listed, not scored.** Turning a reference name like `B.vulgatus.1` into
a taxon string to compare against automatically would mean inventing a mapping this
script has no source for. It prints reference name, ASV and each classifier's call side
by side, 25 rows, and a person reads them.

Standard library only. Needs Python 3.8 or later, Linux or WSL, and conda with the
QIIME 2 environment (see workflow/setup_envs.sh).

Classifiers are not shipped with this repository. Get them from the QIIME 2 Library data
resources, or build one for your own amplicon region with RESCRIPt; either way record
which file you used, since the log captures its UUID.

Example:
    python workflow/05_taxonomy.py -r ~/research/baxter2016/q2/dada2/rep_seqs.qza \\
        -b ~/research/baxter2016/q2/dada2/table.qza \\
        -c gg2=~/ref/gg2-2024.09-515f-806r.qza -c silva=~/ref/silva-138.2-v4.qza \\
        -o ~/research/baxter2016/q2/taxonomy
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import platform
import shlex
import signal
import subprocess
import sys
import zipfile
from datetime import datetime, timezone

__version__ = "0.1.0"

# The rank prefixes SILVA, Greengenes2 and GTDB all use in QIIME 2 artifacts.
RANKS = [("d", "domain"), ("p", "phylum"), ("c", "class"), ("o", "order"),
         ("f", "family"), ("g", "genus"), ("s", "species")]
# Older references write k__ for the top rank instead of d__.
TOP_ALIASES = {"k": "d"}
EMPTY = {"", "unassigned", "unclassified", "uncultured", "metagenome"}

log = logging.getLogger("taxonomy")


class TaxonomyError(RuntimeError):
    """A check failed. The message says which one and why."""


def run_cmd(cmd: list[str], timeout: int) -> tuple[int, str, str]:
    """Run a command with stdin closed and a hard time limit."""
    log.info("run: %s", shlex.join(cmd))
    try:
        proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, start_new_session=True)
    except FileNotFoundError as exc:
        raise TaxonomyError(f"cannot start {cmd[0]}: {exc}") from exc
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.communicate()
        raise TaxonomyError(f"timed out after {timeout} s and was killed: {shlex.join(cmd)}")
    return proc.returncode, out, err


def _tail(text: str, n: int = 15) -> str:
    return "\n".join(text.strip().splitlines()[-n:])


# ---- artifacts -----------------------------------------------------------

def artifact_meta(path: str) -> dict[str, str]:
    """Read type and uuid out of an artifact's metadata.yaml, without QIIME 2."""
    try:
        with zipfile.ZipFile(path) as zf:
            hits = [n for n in zf.namelist()
                    if n.endswith("/metadata.yaml") and n.count("/") == 2]
            if not hits:
                hits = [n for n in zf.namelist() if n.endswith("/metadata.yaml")]
            if not hits:
                raise TaxonomyError(f"{path}: no metadata.yaml, is this a QIIME 2 artifact?")
            text = zf.read(sorted(hits)[0]).decode("utf-8")
    except zipfile.BadZipFile as exc:
        raise TaxonomyError(f"{path} is not a readable artifact: {exc}") from exc
    meta = {}
    for line in text.splitlines():
        if ":" in line:
            key, _, value = line.partition(":")
            meta[key.strip()] = value.strip().strip("'\"")
    return meta


def check_classifier(name: str, path: str) -> dict[str, str]:
    """A classifier must be a TaxonomicClassifier. Anything else wastes hours."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise TaxonomyError(f"classifier {name}: not found or empty: {path}")
    meta = artifact_meta(path)
    kind = meta.get("type", "")
    if "TaxonomicClassifier" not in kind:
        raise TaxonomyError(
            f"classifier {name} is a {kind or 'unknown type'}, not a TaxonomicClassifier. "
            f"Check {path}")
    return meta


def parse_classifiers(values: list[str]) -> dict[str, str]:
    """Read the NAME=PATH pairs, so every output carries a name a human chose."""
    out: dict[str, str] = {}
    for value in values:
        name, sep, path = value.partition("=")
        if not sep or not name.strip() or not path.strip():
            raise TaxonomyError(f"--classifier wants NAME=PATH, got {value!r}")
        name = name.strip()
        if name in out:
            raise TaxonomyError(f"classifier name {name!r} given more than once")
        out[name] = os.path.abspath(os.path.expanduser(path.strip()))
    if not out:
        raise TaxonomyError("no classifiers given")
    return out


# ---- taxonomy strings ----------------------------------------------------

def split_taxon(taxon: str) -> dict[str, str]:
    """Split a taxon string into ranks. Prefixed labels win; otherwise use position."""
    parts = [p.strip() for p in (taxon or "").split(";") if p.strip()]
    by_rank: dict[str, str] = {}
    prefixed = 0
    for i, part in enumerate(parts):
        if len(part) > 3 and part[1:3] == "__":
            code = TOP_ALIASES.get(part[0], part[0])
            label = part[3:].strip()
            prefixed += 1
        elif len(part) == 3 and part[1:3] == "__":
            code, label, prefixed = TOP_ALIASES.get(part[0], part[0]), "", prefixed + 1
        else:
            code = RANKS[i][0] if i < len(RANKS) else None
            label = part
        if code is None:
            continue
        name = dict(RANKS).get(code)
        if name and label.lower() not in EMPTY:
            by_rank[name] = label
    return by_rank


def parse_taxonomy(text: str, name: str) -> dict[str, dict]:
    """Read a QIIME 2 taxonomy.tsv into {feature: {taxon, confidence, ranks}}."""
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise TaxonomyError(f"{name}: no rows")
    header = list(rows[0].keys())
    if header[0] != "Feature ID":
        raise TaxonomyError(f"{name}: first column is {header[0]!r}, expected 'Feature ID'")
    if "Taxon" not in header:
        raise TaxonomyError(f"{name}: no Taxon column")
    out: dict[str, dict] = {}
    for row in rows:
        fid = row["Feature ID"]
        if str(fid).startswith("#q2:types"):
            continue
        if fid in out:
            raise TaxonomyError(f"{name}: feature {fid} appears more than once")
        conf = row.get("Confidence")
        try:
            conf = float(conf) if conf not in (None, "") else None
        except ValueError:
            conf = None
        out[fid] = {"taxon": row["Taxon"], "confidence": conf,
                    "ranks": split_taxon(row["Taxon"])}
    if not out:
        raise TaxonomyError(f"{name}: the table holds no features")
    return out


def read_rep_seq_ids(path: str) -> list[str]:
    """Feature ids from an exported rep_seqs FASTA."""
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    if not ids:
        raise TaxonomyError(f"{path}: no sequences")
    return ids


# ---- reporting -----------------------------------------------------------

def coverage(tax: dict[str, dict], reads: dict[str, float] | None) -> dict[str, dict]:
    """How far down each rank the classifier got, by ASV and by read."""
    total_asvs = len(tax)
    total_reads = sum(reads.values()) if reads else 0.0
    out = {}
    for _, rank in RANKS:
        hit = [f for f, t in tax.items() if rank in t["ranks"]]
        hit_reads = sum(reads.get(f, 0.0) for f in hit) if reads else 0.0
        out[rank] = {
            "asvs": len(hit),
            "asv_fraction": len(hit) / total_asvs if total_asvs else 0.0,
            "reads": hit_reads,
            "read_fraction": hit_reads / total_reads if total_reads else 0.0,
        }
    return out


def label_differences(a: dict[str, dict], b: dict[str, dict]) -> dict[str, dict]:
    """Per rank: of the features both name, how often is the label not the same string.

    This is a string comparison and nothing more. It does NOT measure how often two
    references place an organism differently, and it must not be reported as if it did.

    Greengenes2 writes GTDB names, so the phylum SILVA calls Firmicutes it calls
    Bacillota_A_368345. Those are the same clade under two nomenclatures, and they count
    here as a different label. On a real GG2-vs-SILVA run that put the figure at 85% at
    phylum, essentially all of it naming rather than placement.

    Resolving nomenclature would need a curated synonym table mapping every pair of
    reference vocabularies, kept current as they change. That is deliberately not done:
    a partial table would silently miscount every pair it does not know, and claiming
    agreement we cannot establish is worse than reporting a number that says what it is.

    The fix at the source is to use one reference, which is this stage's default.
    """
    out = {}
    for _, rank in RANKS:
        both = [f for f in a if f in b and rank in a[f]["ranks"] and rank in b[f]["ranks"]]
        differ = [f for f in both if a[f]["ranks"][rank] != b[f]["ranks"][rank]]
        only_a = sum(1 for f in a if rank in a[f]["ranks"] and
                     (f not in b or rank not in b[f]["ranks"]))
        only_b = sum(1 for f in b if rank in b[f]["ranks"] and
                     (f not in a or rank not in a[f]["ranks"]))
        out[rank] = {"compared": len(both), "different_label": len(differ),
                     "different_label_fraction": len(differ) / len(both) if both else 0.0,
                     "only_first": only_a, "only_second": only_b,
                     "examples": differ[:5]}
    return out


def write_coverage(path: str, per_classifier: dict[str, dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["classifier", "rank", "asvs", "asv_fraction", "reads", "read_fraction"])
        for name, cov in per_classifier.items():
            for _, rank in RANKS:
                c = cov[rank]
                w.writerow([name, rank, c["asvs"], round(c["asv_fraction"], 6),
                            round(c["reads"], 1), round(c["read_fraction"], 6)])


def write_label_differences(path: str, pairs: dict[tuple[str, str], dict]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["first", "second", "rank", "compared", "different_label",
                    "different_label_fraction",
                    "named_only_by_first", "named_only_by_second"])
        for (a, b), per_rank in pairs.items():
            for _, rank in RANKS:
                d = per_rank[rank]
                w.writerow([a, b, rank, d["compared"], d["different_label"],
                            round(d["different_label_fraction"], 6),
                            d["only_first"], d["only_second"]])


def write_calls(path: str, features: list[str], tax: dict[str, dict[str, dict]],
                extra: dict[str, str] | None = None) -> None:
    """One row per feature, one column per classifier. Used for the full set and mocks."""
    names = list(tax)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        head = ["feature"] + (["reference_names"] if extra else [])
        for name in names:
            head += [f"{name}_taxon", f"{name}_confidence"]
        w.writerow(head)
        for fid in features:
            row = [fid] + ([extra.get(fid, "")] if extra else [])
            for name in names:
                rec = tax[name].get(fid)
                row += ["" if rec is None else rec["taxon"],
                        "" if rec is None or rec["confidence"] is None
                        else round(rec["confidence"], 4)]
            w.writerow(row)


# ---- exports -------------------------------------------------------------

def export(src: str, dest: str, env: str, timeout: int, expect: str) -> str:
    rc, _, err = run_cmd(["conda", "run", "-n", env, "qiime", "tools", "export",
                          "--input-path", src, "--output-path", dest], timeout)
    if rc != 0:
        raise TaxonomyError(f"qiime tools export failed on {src} with code {rc}:\n{_tail(err)}")
    path = os.path.join(dest, expect)
    if not os.path.isfile(path):
        raise TaxonomyError(f"export finished but wrote no {path}")
    return path


def read_reads_per_feature(table_qza: str, workdir: str, env: str,
                           timeout: int) -> dict[str, float]:
    """Total reads per feature, so coverage can be weighted by abundance."""
    biom = export(table_qza, os.path.join(workdir, "table"), env, timeout,
                  "feature-table.biom")
    tsv = os.path.join(os.path.dirname(biom), "feature-table.tsv")
    rc, _, err = run_cmd(["conda", "run", "-n", env, "biom", "convert", "-i", biom,
                          "-o", tsv, "--to-tsv"], timeout)
    if rc != 0:
        raise TaxonomyError(f"biom convert exited with code {rc}:\n{_tail(err)}")
    totals: dict[str, float] = {}
    with open(tsv) as fh:
        for i, line in enumerate(fh):
            if line.startswith("# Const") or not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            if cells[0].startswith("#OTU"):
                continue
            try:
                totals[cells[0]] = sum(float(c) for c in cells[1:])
            except ValueError as exc:
                raise TaxonomyError(f"{tsv}: {cells[0]} has a non-numeric count") from exc
    if not totals:
        raise TaxonomyError(f"{tsv}: no features")
    return totals


def read_mock_targets(path: str) -> dict[str, str]:
    """Map target sequence -> reference names, from 04_mock.py's targets.tsv."""
    out = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq = (row.get("sequence") or "").strip()
            if seq and not seq.startswith("NO REGION"):
                out[seq] = row.get("reference_names", "")
    if not out:
        raise TaxonomyError(f"{path}: no target sequences. Is this 04_mock's targets.tsv?")
    return out


def read_fasta(path: str) -> dict[str, str]:
    records, name, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if name is not None:
        records[name] = "".join(parts).upper()
    return records


def setup_logging(outdir: str | None) -> None:
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(fmt)
    log.addHandler(console)
    if outdir:
        fh = logging.FileHandler(os.path.join(outdir, "taxonomy_log.txt"), mode="a")
        fh.setFormatter(fmt)
        log.addHandler(fh)


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Classify ASVs against a reference taxonomy and report how deep "
                    "the names go. A second classifier is a diagnostic, not a default.")
    p.add_argument("-r", "--rep-seqs", required=True, help="ASV sequences, e.g. rep_seqs.qza")
    p.add_argument("-c", "--classifier", action="append", default=[], metavar="NAME=PATH",
                   help="a trained classifier and the name to report it under. "
                        "Repeat for more than one")
    p.add_argument("-o", "--outdir", required=True, help="output directory")
    p.add_argument("-b", "--table", help="feature table, to weight coverage by reads")
    p.add_argument("--mock-targets", help="targets.tsv from 04_mock.py, to list what each "
                                          "classifier called the mock's exact matches")
    p.add_argument("--n-jobs", type=int, default=1,
                   help="classify-sklearn jobs. QIIME 2's own default is 1, so this is "
                        "always passed explicitly (default 1)")
    p.add_argument("--confidence", type=float, default=0.7,
                   help="classify-sklearn confidence, default 0.7, which is also QIIME 2's")
    p.add_argument("--env", default="qiime2-amplicon-2025.7", help="conda env with QIIME 2")
    p.add_argument("--timeout", type=int, default=86400,
                   help="seconds before a classify step is killed (default 86400)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p.parse_args(argv)


def run(args) -> None:
    rep_seqs = os.path.abspath(os.path.expanduser(args.rep_seqs))
    outdir = os.path.abspath(os.path.expanduser(args.outdir))
    if not os.path.isfile(rep_seqs) or os.path.getsize(rep_seqs) == 0:
        raise TaxonomyError(f"rep-seqs not found or empty: {rep_seqs}")
    if args.n_jobs == 0 or args.n_jobs < -1:
        raise TaxonomyError(f"--n-jobs must be -1 or 1 or more, got {args.n_jobs}")
    if not 0 <= args.confidence <= 1:
        raise TaxonomyError(f"--confidence must be between 0 and 1, got {args.confidence}")
    classifiers = parse_classifiers(args.classifier)
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir)
    log.info("=" * 70)
    log.info("taxonomy %s | Python %s | %s", __version__, platform.python_version(),
             platform.platform())
    log.info("command: %s", shlex.join([sys.executable] + sys.argv))
    log.info("started: %s", datetime.now(timezone.utc).isoformat(timespec="seconds"))

    for name, path in classifiers.items():
        meta = check_classifier(name, path)
        log.info("classifier %s: %s, uuid %s, %.0f MB", name, meta.get("type", "?"),
                 meta.get("uuid", "?"), os.path.getsize(path) / 1e6)

    fasta = export(rep_seqs, os.path.join(outdir, "seqs"), args.env, args.timeout,
                   "dna-sequences.fasta")
    feature_ids = read_rep_seq_ids(fasta)
    log.info("%d ASVs to classify", len(feature_ids))

    reads = None
    if args.table:
        reads = read_reads_per_feature(os.path.abspath(os.path.expanduser(args.table)),
                                       outdir, args.env, args.timeout)
        missing = [f for f in feature_ids if f not in reads]
        if missing:
            raise TaxonomyError(f"{len(missing)} ASV(s) are in rep-seqs but not the table, "
                                f"for example {missing[0]}. Do they come from the same run?")
        log.info("coverage will also be weighted by reads (%d features in the table)",
                 len(reads))

    tax: dict[str, dict[str, dict]] = {}
    for name, path in classifiers.items():
        qza = os.path.join(outdir, f"taxonomy_{name}.qza")
        rc, _, err = run_cmd(
            ["conda", "run", "-n", args.env, "qiime", "feature-classifier",
             "classify-sklearn", "--i-reads", rep_seqs, "--i-classifier", path,
             "--p-n-jobs", str(args.n_jobs), "--p-confidence", str(args.confidence),
             "--o-classification", qza], args.timeout)
        if rc != 0:
            raise TaxonomyError(f"classify-sklearn failed for {name} with code {rc}:"
                                f"\n{_tail(err)}")
        if not os.path.isfile(qza):
            raise TaxonomyError(f"classify-sklearn finished but wrote no {qza}")
        tsv = export(qza, os.path.join(outdir, f"export_{name}"), args.env, args.timeout,
                     "taxonomy.tsv")
        with open(tsv) as fh:
            tax[name] = parse_taxonomy(fh.read(), tsv)
        unseen = [f for f in feature_ids if f not in tax[name]]
        if unseen:
            raise TaxonomyError(f"{name}: {len(unseen)} ASV(s) got no row at all, "
                                f"for example {unseen[0]}")
        if not any("domain" in t["ranks"] for t in tax[name].values()):
            raise TaxonomyError(f"{name}: not one ASV was placed even at domain level. "
                                "That is a broken classifier or the wrong region, not a "
                                "hard dataset")

    per_classifier = {name: coverage(t, reads) for name, t in tax.items()}
    write_coverage(os.path.join(outdir, "taxonomy_coverage.tsv"), per_classifier)
    write_calls(os.path.join(outdir, "taxonomy_calls.tsv"), feature_ids, tax)

    for name, cov in per_classifier.items():
        parts = []
        for _, rank in RANKS:
            c = cov[rank]
            parts.append(f"{rank} {100 * c['asv_fraction']:.1f}%"
                         + (f"/{100 * c['read_fraction']:.1f}%" if reads else ""))
        log.info("%s classified (by ASV%s): %s", name, "/by read" if reads else "",
                 ", ".join(parts))

    pairs = {}
    names = list(tax)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pairs[(a, b)] = label_differences(tax[a], tax[b])
    if pairs:
        write_label_differences(os.path.join(outdir, "taxonomy_label_differences.tsv"), pairs)
        log.info("more than one classifier given, so label differences are reported below. "
                 "READ THESE AS STRING COMPARISONS, NOT AS PLACEMENT CONFLICTS: two "
                 "references using different nomenclature (GTDB Bacillota_A_368345 vs "
                 "Firmicutes) count as different here while naming the same clade")
        for (a, b), per_rank in pairs.items():
            for _, rank in RANKS:
                d = per_rank[rank]
                if d["compared"]:
                    log.info("%s vs %s at %s: %d compared, %d with a different label "
                             "(%.1f%%), named only by %s %d, only by %s %d", a, b, rank,
                             d["compared"], d["different_label"],
                             100 * d["different_label_fraction"],
                             a, d["only_first"], b, d["only_second"])
    else:
        log.info("one classifier, which is the default and the right one for an analysis. "
                 "The coverage table is the result to report: how deep the names go, by "
                 "ASV and by read. Pass a second classifier only as a diagnostic, to see "
                 "how much of the naming depends on the reference you chose")

    if args.mock_targets:
        targets = read_mock_targets(os.path.abspath(os.path.expanduser(args.mock_targets)))
        seqs = read_fasta(fasta)
        found = {f: targets[s] for f, s in seqs.items() if s in targets}
        log.info("mock: %d of %d reference targets are present as an exact ASV",
                 len(found), len(targets))
        write_calls(os.path.join(outdir, "taxonomy_mock_calls.tsv"),
                    sorted(found), tax, extra=found)
        log.info("mock calls listed in taxonomy_mock_calls.tsv. They are listed, not "
                 "scored: turning a name like B.vulgatus.1 into a taxon string to compare "
                 "against would mean inventing a mapping this script has no source for")

    log.info("note: 253 bp of V4 does not carry species-level information for many taxa "
             "(S. aureus and S. epidermidis are identical over this region in the mock "
             "reference), so read a low species rate as a property of the amplicon")
    log.info("done: %s", os.path.join(outdir, "taxonomy_coverage.tsv"))


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        run(args)
    except TaxonomyError as exc:
        if not log.handlers:
            setup_logging(None)
        log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
