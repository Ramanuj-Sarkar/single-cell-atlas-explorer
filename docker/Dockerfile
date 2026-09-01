# Container image for the PBMC3k Scanpy pipeline (used by `-profile docker`).
# The analysis scripts are baked in at /opt/pipeline so they are available
# inside the container regardless of the host checkout location.

FROM python:3.12-slim

LABEL org.opencontainers.image.title="scanpy-pbmc3k" \
      org.opencontainers.image.description="Scanpy single-cell RNA-seq analysis environment for the 10x PBMC3k dataset (QC, normalization, dimensionality reduction, clustering, annotation, MultiQC)" \
      org.opencontainers.image.version="1.0.0" \
      org.opencontainers.image.licenses="MIT"

# OpenMP runtime (needed by some numerical binaries); build tools are not
# required because all dependencies ship manylinux wheels for cp312.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgomp1 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY docker/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r /tmp/requirements.txt

COPY obtain_dataset.py /opt/pipeline/obtain_dataset.py
COPY analysis/ /opt/pipeline/analysis/
# COPY preserves host permissions; ensure the (possibly 600-mode) scripts are
# readable by the non-root runtime user.
RUN chmod -R a+rX /opt/pipeline

WORKDIR /opt/pipeline

# keep caches inside the container (writable, and avoids stale host caches)
ENV NUMBA_CACHE_DIR=/tmp/numba-cache \
    MPLCONFIGDIR=/tmp/matplotlib-cache \
    PYTHONUNBUFFERED=1

# run as a non-root user
RUN useradd --create-home pipeline
USER pipeline
