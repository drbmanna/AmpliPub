"""Tests for scripts/check_env_lock.py and the committed lockfiles.

The check is only worth having if it fires, so each kind of difference (a version change,
a rebuild with the same version, a missing package, an extra package, another platform)
is shown to be caught. conda itself is replaced by a fake that prints a given listing.
"""

from __future__ import annotations

import csv
import importlib.util
import pathlib
import sys

import pytest

HERE = pathlib.Path(__file__).resolve().parent
ENVS = HERE.parent / "envs"
spec = importlib.util.spec_from_file_location("check_env_lock",
                                              HERE.parent / "scripts" / "check_env_lock.py")
cel = importlib.util.module_from_spec(spec)
sys.modules["check_env_lock"] = cel
spec.loader.exec_module(cel)

CF = "https://conda.anaconda.org/conda-forge/linux-64"
MD5A, MD5B, MD5C = "a" * 32, "b" * 32, "c" * 32


def listing(*lines: str, platform: str = "linux-64") -> str:
    return "\n".join(["# This file may be used to create an environment using:",
                      f"# platform: {platform}", "@EXPLICIT", *lines]) + "\n"


LOCK = listing(f"{CF}/r-rbiom-2.2.1-r45h3697838_1.conda#{MD5A}",
               f"{CF}/r-base-4.5.3-h1_0.conda#{MD5B}",
               f"{CF}/python-dateutil-2.9.0.post0-pyhe01879c_2.conda#{MD5C}")


def run(tmp_path, monkeypatch, installed: str, *extra: str) -> tuple[int, list[dict]]:
    lock = tmp_path / "env.lock"
    lock.write_text(LOCK)
    monkeypatch.setattr(cel, "conda_explicit", lambda env: installed)
    out = tmp_path / "diff.tsv"
    code = cel.main(["--env", "x", "--lock", str(lock), "-o", str(out), *extra])
    with open(out) as fh:
        return code, list(csv.DictReader(fh, delimiter="\t"))


def test_identical_env_passes_with_an_empty_table(tmp_path, monkeypatch):
    code, rows = run(tmp_path, monkeypatch, LOCK)
    assert code == 0 and rows == []


def test_order_of_lines_does_not_matter(tmp_path, monkeypatch):
    head, body = LOCK.split("@EXPLICIT\n")
    shuffled = head + "@EXPLICIT\n" + "\n".join(reversed(body.strip().splitlines())) + "\n"
    assert run(tmp_path, monkeypatch, shuffled)[0] == 0


def test_the_rbiom_case_a_different_version_fails(tmp_path, monkeypatch):
    bad = LOCK.replace("r-rbiom-2.2.1-r45h3697838_1", "r-rbiom-3.1.0-r45h0_0")
    code, rows = run(tmp_path, monkeypatch, bad)
    assert code == 1
    assert rows == [{"package": "r-rbiom", "status": "different",
                     "locked": "r-rbiom-2.2.1-r45h3697838_1.conda",
                     "installed": "r-rbiom-3.1.0-r45h0_0.conda"}]


def test_same_file_name_with_another_md5_fails(tmp_path, monkeypatch):
    code, rows = run(tmp_path, monkeypatch, LOCK.replace(f"#{MD5A}", "#" + "d" * 32))
    assert code == 1 and rows[0]["package"] == "r-rbiom"


def test_missing_and_extra_packages_are_both_reported(tmp_path, monkeypatch):
    inst = LOCK.replace(f"{CF}/r-base-4.5.3-h1_0.conda#{MD5B}\n", "")
    inst += f"{CF}/r-vegan-2.7_5-r45h1_0.conda#{MD5B}\n"
    code, rows = run(tmp_path, monkeypatch, inst)
    assert code == 1
    assert {(r["package"], r["status"]) for r in rows} == {
        ("r-base", "missing from env"), ("r-vegan", "not in lockfile")}


def test_another_platform_fails(tmp_path, monkeypatch):
    code, rows = run(tmp_path, monkeypatch, LOCK.replace("linux-64", "osx-arm64", 1))
    assert code == 1 and rows[0]["package"] == "(platform)"


def test_report_only_records_the_difference_but_exits_zero(tmp_path, monkeypatch, capsys):
    bad = LOCK.replace("r-rbiom-2.2.1-r45h3697838_1", "r-rbiom-3.1.0-r45h0_0")
    code, rows = run(tmp_path, monkeypatch, bad, "--report-only")
    assert code == 0 and len(rows) == 1
    assert "WARNING" in capsys.readouterr().err


def test_names_with_hyphens_are_split_correctly():
    assert cel.package_key("python-dateutil-2.9.0.post0-pyhe01879c_2.conda") == "python-dateutil"
    assert cel.package_key("libgcc-ng-15.1.0-h69a702a_4.tar.bz2") == "libgcc-ng"


@pytest.mark.parametrize("text, message", [
    ("# platform: linux-64\nhttps://x/a-1-0.conda#" + MD5A + "\n", "no @EXPLICIT"),
    (listing("https://x/a-1-0.conda"), "no md5"),
    (listing(), "no packages"),
    (listing(f"https://x/a-1-0.conda#{MD5A}", f"https://y/a-1-0.conda#{MD5B}"), "twice"),
])
def test_malformed_listings_are_refused(text, message):
    with pytest.raises(cel.LockError, match=message):
        cel.parse_explicit(text, "t")


def test_unreadable_env_exits_2(tmp_path, monkeypatch):
    lock = tmp_path / "env.lock"
    lock.write_text(LOCK)

    def boom(env):
        raise cel.LockError("conda list -n x failed")
    monkeypatch.setattr(cel, "conda_explicit", boom)
    assert cel.main(["--env", "x", "--lock", str(lock), "-o", str(tmp_path / "d.tsv")]) == 2


@pytest.mark.parametrize("name", ["qiime2-amplicon-2025.7", "amplipub-qc",
                                  "amplipub-snakemake", "amplipub-r"])
def test_committed_lockfiles_are_complete(name):
    pkgs = cel.parse_explicit((ENVS / f"{name}.lock").read_text(), name)
    assert pkgs["__platform__"] == "linux-64"
    assert len(pkgs) > 50


def test_r_lockfile_holds_the_rbiom_pin():
    pkgs = cel.parse_explicit((ENVS / "amplipub-r.lock").read_text(), "r")
    rbiom = [f for f in pkgs if cel.package_key(f) == "r-rbiom"]
    assert len(rbiom) == 1 and rbiom[0].startswith("r-rbiom-2.2.1-")
