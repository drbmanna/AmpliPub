#!/usr/bin/env bash
# Create the conda envs the workflow stages call, from the committed lockfiles.
#
#   qiime2-amplicon-2025.7   QIIME 2 amplicon release       envs/qiime2-amplicon-2025.7.lock
#   amplipub-qc              FastQC + MultiQC               envs/amplipub-qc.lock
#   amplipub-snakemake       Snakemake                      envs/amplipub-snakemake.lock
#   amplipub-r               R and every package AmpliPub   envs/amplipub-r.lock
#                            uses, with AmpliPub itself installed from this checkout
#
# The lockfiles (`conda list --explicit --md5`, linux-64) pin every package, dependencies
# included, so no solver runs and the env is the one the tests passed on. The *.yml files
# pin only the packages we name; on 2026-09-13 the solver filled the gap with an rbiom that
# mia could not load. The yml files are the readable spec, used only to upgrade on purpose.
#
# Usage:
#   bash workflow/setup_envs.sh              create missing envs from the lockfiles, and
#                                            refuse if an existing env differs from its lock
#   bash workflow/setup_envs.sh --rebuild    remove and recreate every env from its lock
#   bash workflow/setup_envs.sh --from-spec  create missing envs by solving the yml files
#                                            (deliberate upgrade); re-lock afterwards with
#                                            conda list -n ENV --explicit --md5 > envs/ENV.lock
set -euo pipefail

mode=lock
case "${1:-}" in
  "") ;;
  --rebuild) mode=rebuild ;;
  --from-spec) mode=spec ;;
  *) echo "ERROR: unknown option $1 (use --rebuild or --from-spec)" >&2; exit 2 ;;
esac

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/.." && pwd)"
q2_env="qiime2-amplicon-2025.7"
qc_env="amplipub-qc"
smk_env="amplipub-snakemake"
r_env="amplipub-r"
# name|spec|lock. The QIIME 2 spec is the official release file.
envs=(
  "$q2_env|https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml|$here/envs/$q2_env.lock"
  "$qc_env|$here/envs/qc.yml|$here/envs/$qc_env.lock"
  "$smk_env|$here/envs/snakemake.yml|$here/envs/$smk_env.lock"
  "$r_env|$here/envs/amplipub-r.yml|$here/envs/$r_env.lock"
)

command -v conda >/dev/null || { echo "ERROR: conda not found on PATH" >&2; exit 1; }
if [ "$mode" != spec ] && [ "$(uname -s)-$(uname -m)" != "Linux-x86_64" ]; then
  echo "ERROR: the lockfiles are linux-64; on this platform use --from-spec" >&2; exit 1
fi
for e in "${envs[@]}"; do
  IFS='|' read -r _ spec lock <<<"$e"
  if [ "$mode" = spec ]; then need="$spec"; else need="$lock"; fi
  if [[ "$need" != http* ]] && [ ! -s "$need" ]; then
    echo "ERROR: missing $need" >&2; exit 1
  fi
done

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }
check_py="$here/scripts/check_env_lock.py"
drift_dir="$(mktemp -d)"
trap 'rm -rf "$drift_dir"' EXIT

# Env creation can be slow but must never hang forever; stdin closed so nothing waits on input.
for e in "${envs[@]}"; do
  IFS='|' read -r name spec lock <<<"$e"
  if env_exists "$name" && [ "$mode" = rebuild ]; then
    echo "removing: $name"
    timeout 1800 conda env remove -y -n "$name" </dev/null
  fi
  if env_exists "$name"; then
    echo "exists: $name"
    if [ "$mode" != spec ]; then
      python3 "$check_py" --env "$name" --lock "$lock" -o "$drift_dir/$name.tsv" || {
        echo "ERROR: $name is not the locked env. Fix with --rebuild." >&2; exit 1; }
    fi
  elif [ "$mode" = spec ]; then
    echo "creating from spec: $name"
    timeout 7200 conda env create -n "$name" --file "$spec" </dev/null
    echo "NOTE: $name was solved, not locked. Test it, then re-lock and commit envs/$name.lock"
  else
    echo "creating from lock: $name"
    timeout 7200 conda create -y -n "$name" --file "$lock" </dev/null
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
