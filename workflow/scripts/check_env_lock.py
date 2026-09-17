#!/usr/bin/env python3
"""Compare an installed conda env with its committed lockfile, package by package.

The env files in workflow/envs/*.yml pin only the packages AmpliPub names. Everything
those pull in is left to the solver, and on 2026-09-13 the solver picked an rbiom that
mia 1.18.0 could not load. The lockfiles (`conda list --explicit --md5`) pin every
package, down to the build and its md5. This script says whether the env a run actually
used is the one in the lockfile, and writes every difference to a table.

A package is identified by its full URL line, md5 included, so a rebuilt package with
the same version still counts as a difference.

Exit status: 0 when the env matches the lockfile, 1 when it does not, 2 when the check
could not be made. `--report-only` records a mismatch without failing on it.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys

__version__ = "0.1.0"

CONDA_TIMEOUT = 600


class LockError(Exception):
    pass


def parse_explicit(text: str, name: str) -> dict[str, str]:
    """Package file name -> full URL line, from `conda list --explicit --md5` output."""
    lines = [ln.strip() for ln in text.splitlines()]
    if "@EXPLICIT" not in lines:
        raise LockError(f"{name}: no @EXPLICIT line, so this is not an explicit lockfile")
    platforms = [ln.split(":", 1)[1].strip() for ln in lines if ln.startswith("# platform:")]
    pkgs: dict[str, str] = {}
    for ln in lines:
        if not ln or ln.startswith("#") or ln == "@EXPLICIT":
            continue
        if "#" not in ln:
            raise LockError(f"{name}: {ln!r} has no md5, export with --md5")
        fname = ln.split("#", 1)[0].rsplit("/", 1)[-1]
        if fname in pkgs:
            raise LockError(f"{name}: {fname} is listed twice")
        pkgs[fname] = ln
    if not pkgs:
        raise LockError(f"{name}: lists no packages")
    return {"__platform__": platforms[0] if platforms else "", **pkgs}


def package_key(fname: str) -> str:
    """Package name without version and build: 'r-rbiom-2.2.1-r45h_1.conda' -> 'r-rbiom'."""
    stem = fname.removesuffix(".conda").removesuffix(".tar.bz2")
    return stem.rsplit("-", 2)[0]


def compare(locked: dict[str, str], installed: dict[str, str]) -> list[dict[str, str]]:
    """One row per package that is missing, extra, or different."""
    lock_by_name = {package_key(f): f for f in locked if f != "__platform__"}
    inst_by_name = {package_key(f): f for f in installed if f != "__platform__"}
    rows = []
    if locked["__platform__"] != installed["__platform__"]:
        rows.append({"package": "(platform)", "status": "different",
                     "locked": locked["__platform__"], "installed": installed["__platform__"]})
    for name in sorted(set(lock_by_name) | set(inst_by_name)):
        lf, inf = lock_by_name.get(name), inst_by_name.get(name)
        if lf is None:
            rows.append({"package": name, "status": "not in lockfile",
                         "locked": "", "installed": inf})
        elif inf is None:
            rows.append({"package": name, "status": "missing from env",
                         "locked": lf, "installed": ""})
        elif locked[lf] != installed[inf]:
            rows.append({"package": name, "status": "different",
                         "locked": lf, "installed": inf})
    return rows


def conda_explicit(env: str) -> str:
    cmd = ["conda", "list", "-n", env, "--explicit", "--md5"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=CONDA_TIMEOUT,
                              stdin=subprocess.DEVNULL)
    except FileNotFoundError as exc:
        raise LockError("conda is not on PATH") from exc
    except subprocess.TimeoutExpired as exc:
        raise LockError(f"conda list timed out after {CONDA_TIMEOUT} s") from exc
    if proc.returncode != 0:
        raise LockError(f"conda list -n {env} failed: {proc.stderr.strip()[-500:]}")
    return proc.stdout


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--env", required=True, help="name of the installed conda env")
    p.add_argument("--lock", required=True, help="the committed explicit lockfile")
    p.add_argument("-o", "--output", required=True, help="TSV of differences (header only if none)")
    p.add_argument("--report-only", action="store_true",
                   help="write the differences but exit 0 even if there are some")
    args = p.parse_args(argv)

    try:
        with open(args.lock) as fh:
            locked = parse_explicit(fh.read(), args.lock)
        installed = parse_explicit(conda_explicit(args.env), f"env {args.env}")
    except (OSError, LockError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    rows = compare(locked, installed)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["package", "status", "locked", "installed"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    n_locked = len(locked) - 1
    if not rows:
        print(f"env {args.env} matches {args.lock}: {n_locked} packages")
        return 0
    print(f"{'WARNING' if args.report_only else 'ERROR'}: env {args.env} differs from "
          f"{args.lock} in {len(rows)} of {n_locked} packages; see {args.output}",
          file=sys.stderr)
    for r in rows[:10]:
        print(f"  {r['package']}: {r['status']} (locked {r['locked'] or '-'}, "
              f"installed {r['installed'] or '-'})", file=sys.stderr)
    if not args.report_only:
        print("Rebuild the env from the lockfile (bash workflow/setup_envs.sh --rebuild), or, "
              "for a deliberate upgrade, re-lock it and commit the new lockfile.",
              file=sys.stderr)
    return 0 if args.report_only else 1


if __name__ == "__main__":
    sys.exit(main())
