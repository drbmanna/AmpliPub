#!/usr/bin/env bash
# Create the conda envs the workflow stages call. Safe to re-run: existing envs are
# left untouched and only checked.
#
#   qiime2-amplicon-2025.7   QIIME 2 amplicon release, from the official env file
#   amplipub-qc              FastQC + MultiQC, from envs/qc.yml
#
# Usage: bash workflow/setup_envs.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
q2_env="qiime2-amplicon-2025.7"
q2_yml="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"
qc_env="amplipub-qc"
qc_yml="$here/envs/qc.yml"

command -v conda >/dev/null || { echo "ERROR: conda not found on PATH" >&2; exit 1; }
[ -s "$qc_yml" ] || { echo "ERROR: missing $qc_yml" >&2; exit 1; }

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

# Env solves can be slow but must never hang forever; stdin closed so nothing waits on input.
if env_exists "$q2_env"; then
  echo "exists: $q2_env"
else
  echo "creating: $q2_env"
  timeout 3600 conda env create -n "$q2_env" --file "$q2_yml" </dev/null
fi

if env_exists "$qc_env"; then
  echo "exists: $qc_env"
else
  echo "creating: $qc_env"
  timeout 3600 conda env create -n "$qc_env" --file "$qc_yml" </dev/null
fi

echo "--- versions"
timeout 120 conda run -n "$q2_env" qiime --version </dev/null | head -1
timeout 120 conda run -n "$qc_env" fastqc --version </dev/null
timeout 120 conda run -n "$qc_env" multiqc --version </dev/null
