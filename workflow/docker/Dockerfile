# AmpliPub workflow image: the four locked conda envs, nothing else.
#
# Base: the QIIME 2 amplicon 2025.7 image, used for its conda. Its own QIIME 2 env has
# the same package versions as our lockfile but 68 of 769 packages are different builds
# (checked 2026-09-17), so the env is recreated from the lockfile rather than trusted.
# The code is not in the image. It is mounted at run time, so the image never goes stale
# and the recorded commit is the checkout that ran.
#
# Build from the repository root:
#   docker build -f workflow/docker/Dockerfile -t amplipub-workflow:2025.7-lock1 workflow/envs
#
# Run as the calling user. The AmpliPub R library is left writable because rule
# install_amplipub installs the mounted checkout into it on every run; the change lives
# only as long as the container.

FROM quay.io/qiime2/amplicon:2025.7

COPY qiime2-amplicon-2025.7.lock amplipub-qc.lock amplipub-snakemake.lock amplipub-r.lock /opt/locks/

RUN set -eux; \
    conda env remove -y -n qiime2-amplicon-2025.7; \
    for e in qiime2-amplicon-2025.7 amplipub-qc amplipub-snakemake amplipub-r; do \
      conda create -y -q -n "$e" --file "/opt/locks/$e.lock"; \
    done; \
    conda clean -y --all; \
    chmod -R a+rwX /opt/conda/envs/amplipub-r/lib/R/library

# The workflow calls `conda run -n <env>` and `snakemake`; put Snakemake's env first.
ENV PATH=/opt/conda/envs/amplipub-snakemake/bin:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    CONDA_PREFIX= \
    HOME=/tmp \
    MPLBACKEND=agg \
    LC_ALL=C.UTF-8 \
    LANG=C.UTF-8

LABEL org.opencontainers.image.title="amplipub-workflow" \
      org.opencontainers.image.description="Locked conda envs for the AmpliPub Snakemake workflow"

CMD ["/bin/bash"]
