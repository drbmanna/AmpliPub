#!/usr/bin/env python3
"""Fetch an input from a URL, or take it from a local path, and verify its sha256.

Every file the workflow did not produce itself enters through here, so each one is
checked against a checksum written in the config before any stage uses it. A download
is written to a temporary name and only renamed once the checksum matches, so a
truncated or substituted file never appears under the final name.

Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
import urllib.request

__version__ = "0.1.0"


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--url", help="download from this URL")
    source.add_argument("--path", help="use this local file")
    p.add_argument("--sha256", required=True, help="expected sha256 of the file")
    p.add_argument("-o", "--output", required=True, help="where the verified file goes")
    p.add_argument("--link", action="store_true",
                   help="for --path: symlink instead of copying, e.g. for a large classifier")
    p.add_argument("--timeout", type=int, default=300, help="seconds for a download (default 300)")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = p.parse_args(argv)

    expected = args.sha256.strip().lower()
    out = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    origin = args.url or os.path.abspath(args.path)

    if args.url:
        tmp = out + ".part"
        with urllib.request.urlopen(args.url, timeout=args.timeout) as resp, open(tmp, "wb") as fh:
            shutil.copyfileobj(resp, fh)
        checked = tmp
    else:
        if not os.path.isfile(args.path):
            print(f"ERROR: {args.path} does not exist", file=sys.stderr)
            return 1
        checked = os.path.abspath(args.path)

    got = sha256_of(checked)
    if got != expected:
        if args.url:
            os.remove(checked)
        print(f"ERROR: sha256 mismatch for {origin}\n  expected {expected}\n  got      {got}",
              file=sys.stderr)
        return 1

    if args.url:
        os.replace(checked, out)
    elif args.link:
        if os.path.lexists(out):
            os.remove(out)
        os.symlink(checked, out)
    else:
        shutil.copyfile(checked, out)

    print(f"verified {out}\n  source {origin}\n  sha256 {got}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
