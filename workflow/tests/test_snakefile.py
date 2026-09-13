"""Tests for the Snakefile and its config schema.

Two layers, so that CI without Snakemake still checks something real:

- The schema tests need only jsonschema and PyYAML. They skip if either is missing.
- The dry-run test needs Snakemake. It builds the full job graph from the template config
  with every path replaced by a temporary one, and runs nothing. It skips if Snakemake is
  not installed.
"""

from __future__ import annotations

import copy
import shutil
import subprocess
from pathlib import Path

import pytest

WORKFLOW = Path(__file__).resolve().parents[1]
SCHEMA = WORKFLOW / "config" / "schema.yaml"
TEMPLATE = WORKFLOW / "config" / "config.template.yaml"
SNAKEFILE = WORKFLOW / "Snakefile"

yaml = pytest.importorskip("yaml")
jsonschema = pytest.importorskip("jsonschema")

ZEROS = "0" * 64


def load(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def validate(cfg: dict) -> None:
    # The validator Snakemake's utils.validate() uses, whatever the schema declares.
    jsonschema.Draft202012Validator(load(SCHEMA)).validate(cfg)


def test_template_is_valid_against_the_schema():
    validate(load(TEMPLATE))


def test_the_schema_is_written_in_the_dialect_snakemake_validates_with():
    # Read from snakemake/utils.py (9.26.1): validate() always builds a
    # Draft202012Validator and warns "No validator found" unless $schema equals that
    # validator's own identifier. Two earlier choices here were wrong: draft-06 https
    # (UnknownDialect) and draft-06 http, which only moved the warning.
    schema = load(SCHEMA)
    assert schema["$schema"] == jsonschema.Draft202012Validator.META_SCHEMA["$schema"]
    jsonschema.Draft202012Validator.check_schema(schema)


@pytest.mark.parametrize("breakage, message", [
    (lambda c: c.pop("outdir"), "outdir"),
    (lambda c: c["classifier"].update(sha256="not-a-checksum"), "sha256"),
    (lambda c: c["references"]["metadata"].pop("path"), "url"),
    (lambda c: c["analysis"].update(normalizations=["tss", "vst"]), "vst"),
    (lambda c: c["threads"].update(dada2=0), "minimum"),
])
def test_broken_configs_are_refused(breakage, message):
    cfg = copy.deepcopy(load(TEMPLATE))
    breakage(cfg)
    with pytest.raises(jsonschema.ValidationError) as err:
        validate(cfg)
    assert message in str(err.value)


def tiny_config(tmp: Path) -> dict:
    """The template with every path pointing somewhere real under tmp."""
    cfg = copy.deepcopy(load(TEMPLATE))
    for name in ("metadata.tsv", "classifier.qza", "mock.fasta"):
        (tmp / name).write_text("placeholder\n")
    cfg["outdir"] = str(tmp / "results")
    cfg["references"]["metadata"] = {"path": str(tmp / "metadata.tsv"), "sha256": ZEROS}
    cfg["references"]["mock"] = {"path": str(tmp / "mock.fasta"), "sha256": ZEROS}
    cfg["classifier"].update(path=str(tmp / "classifier.qza"), sha256=ZEROS)
    cfg["mock"]["samples"] = ["mockA"]
    cfg["resolution"]["run"] = True
    cfg["resolution"]["primers"] = ["v4=GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT"]
    cfg["diversity"]["depth"] = 1000
    cfg["analysis"]["group"] = "group"
    return cfg


@pytest.mark.skipif(shutil.which("snakemake") is None, reason="Snakemake is not installed")
def test_dry_run_schedules_every_stage_in_order(tmp_path):
    cfg_path = tmp_path / "config.yaml"
    with open(cfg_path, "w") as fh:
        yaml.safe_dump(tiny_config(tmp_path), fh)

    proc = subprocess.run(
        ["snakemake", "-s", str(SNAKEFILE), "--configfile", str(cfg_path),
         "--directory", str(tmp_path / "work"), "-n", "--quiet", "rules"],
        capture_output=True, text=True, timeout=300, stdin=subprocess.DEVNULL,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    out = proc.stdout + proc.stderr
    for rule in ["reference_metadata", "reference_mock", "classifier", "fetch_sra", "qc_raw",
                 "import_demux", "primers", "quality", "dada2", "mock", "taxonomy",
                 "resolution", "collapse", "filter", "tree", "diversity",
                 "provenance_environments", "amplipub_analysis", "report"]:
        assert rule in out, f"rule {rule} was not scheduled"
