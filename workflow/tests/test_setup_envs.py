"""Tests for setup_envs.sh, in particular that --rebuild cannot delete someone else's env.

Until 2026-09-18 the QIIME environment was named `qiime2-amplicon-2025.7`, which is the
name QIIME 2's own installation docs use. Any existing QIIME 2 user already had one, so a
clean install refused, and the error pointed at `--rebuild`, which would have deleted it.
The names are now prefixed and ownership is recorded in the env prefix. Both halves are
tested here: conda is replaced by a fake so the guard can be shown to fire without
touching a real environment.
"""

from __future__ import annotations

import os
import pathlib
import subprocess

import pytest

HERE = pathlib.Path(__file__).resolve().parent
WF = HERE.parent
SCRIPT = WF / "setup_envs.sh"
MARKER = ".amplipub-env"

pytestmark = pytest.mark.skipif(os.name == "nt", reason="bash script, POSIX only")

FAKE_CONDA = """#!/usr/bin/env bash
# Minimal conda stand-in. Prints the listing in $FAKE_ENV_LIST and records
# destructive calls in $FAKE_LOG instead of performing them.
if [ "$1" = "env" ] && [ "$2" = "list" ]; then
  cat "$FAKE_ENV_LIST"
  exit 0
fi
if [ "$1" = "env" ] && [ "$2" = "remove" ]; then
  echo "remove $*" >> "$FAKE_LOG"
  exit 0
fi
if [ "$1" = "create" ] || { [ "$1" = "env" ] && [ "$2" = "create" ]; }; then
  echo "create $*" >> "$FAKE_LOG"
  exit 0
fi
exit 0
"""


def make_env(tmp_path, name: str, owned: bool) -> pathlib.Path:
    prefix = tmp_path / "envs" / name
    prefix.mkdir(parents=True)
    if owned:
        (prefix / MARKER).write_text("created_by\tAmpliPub\n", encoding="utf-8")
    return prefix


def run_setup(tmp_path, envs: dict[str, bool], *args: str):
    """Run setup_envs.sh with a fake conda reporting `envs` {name: amplipub_owned}."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    conda = bin_dir / "conda"
    conda.write_text(FAKE_CONDA, encoding="utf-8")
    conda.chmod(0o755)

    listing = tmp_path / "envlist.txt"
    lines = ["# conda environments:", "#", "base                     /opt/conda"]
    for name, owned in envs.items():
        lines.append(f"{name}                  {make_env(tmp_path, name, owned)}")
    listing.write_text("\n".join(lines) + "\n", encoding="utf-8")

    log = tmp_path / "conda.log"
    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}{os.pathsep}{env['PATH']}"
    env["FAKE_ENV_LIST"] = str(listing)
    env["FAKE_LOG"] = str(log)

    proc = subprocess.run(["bash", str(SCRIPT), *args], capture_output=True,
                          text=True, env=env, timeout=120)
    actions = log.read_text(encoding="utf-8") if log.exists() else ""
    return proc, actions


def test_rebuild_refuses_to_delete_an_env_amplipub_did_not_create(tmp_path):
    proc, actions = run_setup(tmp_path, {"amplipub-qiime2-2025.7": False}, "--rebuild")
    assert proc.returncode != 0
    assert "did not create it" in proc.stderr
    assert "Refusing to delete" in proc.stderr
    assert "remove" not in actions, "the guard must fire before conda env remove runs"


def test_rebuild_removes_an_env_amplipub_did_create(tmp_path):
    proc, actions = run_setup(tmp_path, {"amplipub-qiime2-2025.7": True}, "--rebuild")
    assert "remove" in actions
    assert "amplipub-qiime2-2025.7" in actions


def test_drift_on_a_foreign_env_does_not_point_the_user_at_rebuild(tmp_path):
    # The old message said "Fix with --rebuild", which was the command that would
    # have deleted their QIIME 2 install.
    proc, _ = run_setup(tmp_path, {"amplipub-qiime2-2025.7": False})
    assert proc.returncode != 0
    combined = proc.stdout + proc.stderr
    assert "Do not rebuild it" in combined
    assert "Point the config at a different env name" in combined


def test_no_env_name_collides_with_the_qiime2_install_docs_name(tmp_path):
    text = SCRIPT.read_text(encoding="utf-8")
    for line in text.splitlines():
        if line.startswith(("q2_env=", "qc_env=", "smk_env=", "r_env=")):
            value = line.split("=", 1)[1].strip().strip('"')
            assert value.startswith("amplipub-"), f"{value} is not namespaced"
            assert value != "qiime2-amplicon-2025.7"


def test_every_managed_env_has_a_matching_lockfile():
    names = []
    for line in SCRIPT.read_text(encoding="utf-8").splitlines():
        if line.startswith(("q2_env=", "qc_env=", "smk_env=", "r_env=")):
            names.append(line.split("=", 1)[1].strip().strip('"'))
    assert names, "no env names found in setup_envs.sh"
    for n in names:
        assert (WF / "envs" / f"{n}.lock").is_file(), f"missing lockfile for {n}"
