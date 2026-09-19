"""The README title must carry the package version from DESCRIPTION.

A version bump that forgets the README would otherwise go unnoticed, so CI checks it.
"""

from __future__ import annotations

import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parents[2]


def description_version(text: str) -> str:
    m = re.search(r"^Version:\s*(\S+)\s*$", text, re.M)
    assert m, "DESCRIPTION has no Version field"
    return m.group(1)


def readme_title_mismatch(readme: str, version: str) -> str | None:
    """Return a reason if the first line is not '# AmpliPub <version>', else None."""
    first = readme.splitlines()[0] if readme else ""
    expected = f"# AmpliPub {version}"
    return None if first == expected else f"README title is {first!r}, expected {expected!r}"


def test_readme_title_matches_description_version():
    version = description_version((REPO / "DESCRIPTION").read_text(encoding="utf-8"))
    readme = (REPO / "README.md").read_text(encoding="utf-8")
    reason = readme_title_mismatch(readme, version)
    assert reason is None, reason


def test_guard_fires_on_stale_or_missing_version():
    assert readme_title_mismatch("# AmpliPub 0.0.1\n", "0.0.2") is not None
    assert readme_title_mismatch("# AmpliPub\n", "0.0.1") is not None
    assert readme_title_mismatch("", "0.0.1") is not None
    assert readme_title_mismatch("# AmpliPub 0.0.1\ntext\n", "0.0.1") is None
