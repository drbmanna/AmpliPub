#!/usr/bin/env bash
# Create the conda envs the workflow stages call. Safe to re-run: existing envs are
# left untouched and only checked.
#
#   qiime2-amplicon-2025.7   QIIME 2 amplicon release, from the official env file
#   amplipub-qc              FastQC + MultiQC, from envs/qc.yml
#   amplipub-snakemake       Snakemake, from envs/snakemake.yml
#   amplipub-r               R and every package AmpliPub uses, from envs/amplipub-r.yml,
#                            with AmpliPub itself installed from this checkout
#
# Usage: bash workflow/setup_envs.sh
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/.." && pwd)"
q2_env="qiime2-amplicon-2025.7"
q2_yml="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"
qc_env="amplipub-qc"
qc_yml="$here/envs/qc.yml"
smk_env="amplipub-snakemake"
smk_yml="$here/envs/snakemake.yml"
r_env="amplipub-r"
r_yml="$here/envs/amplipub-r.yml"

command -v conda >/dev/null || { echo "ERROR: conda not found on PATH" >&2; exit 1; }
for f in "$qc_yml" "$smk_yml" "$r_yml"; do
  [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

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

for pair in "$smk_env:$smk_yml" "$r_env:$r_yml"; do
  name="${pair%%:*}"; yml="${pair#*:}"
  if env_exists "$name"; then
    echo "exists: $name"
  else
    echo "creating: $name"
    timeout 7200 conda env create -n "$name" --file "$yml" </dev/null
  fi
done

# AmpliPub is installed from this checkout, not from a repository, so the package the
# pipeline runs is exactly the commit that holds the Snakefile. Reinstalled on every run
# of this script, since the checkout may have moved on.
echo "installing AmpliPub from $repo into $r_env"
timeout 1800 conda run -n "$r_env" R CMD INSTALL --no-multiarch "$repo" </dev/null

echo "--- versions"
timeout 120 conda run -n "$q2_env" qiime --version </dev/null | head -1
timeout 120 conda run -n "$qc_env" fastqc --version </dev/null
timeout 120 conda run -n "$qc_env" multiqc --version </dev/null
timeout 120 conda run -n "$smk_env" snakemake --version </dev/null
# Load every package, not just read its version: on 2026-09-13 mia reported 1.18.0 installed
# and still failed to load, because the solver had picked an rbiom it was not built against.
timeout 900 conda run -n "$r_env" Rscript -e '
cat(R.version.string, "\n")
pkgs <- c("AmpliPub", "ape", "Biostrings", "biomformat", "boot", "cli", "ggplot2", "Matrix",
          "mia", "permute", "rbiom", "rlang", "S4Vectors", "scales", "SummarizedExperiment",
          "TreeSummarizedExperiment", "vegan", "ALDEx2", "ANCOMBC", "Maaslin2",
          "MicrobiomeStat", "microbiome", "metagenomeSeq", "lme4", "rstatix", "ggrepel",
          "withr", "zip", "testthat", "devtools", "rmarkdown", "knitr", "yaml")
failed <- character(0)
for (p in pkgs) {
  msg <- tryCatch({
    suppressPackageStartupMessages(loadNamespace(p))
    paste("loads", as.character(utils::packageVersion(p)))
  }, error = function(e) { failed <<- c(failed, p); paste("FAILS:", conditionMessage(e)) })
  cat(sprintf("%-26s %s\n", p, msg))
}
if (length(failed)) { cat("NOT LOADING:", paste(failed, collapse = ", "), "\n"); quit(status = 1) }
' </dev/null
