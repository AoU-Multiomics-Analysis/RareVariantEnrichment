FROM python:3.12-slim-bookworm
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        tabix \
        ca-certificates \
        r-base \
        r-cran-data.table \
        r-cran-ggplot2 \
        r-cran-ggrepel \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/rare-variant-enrichment
COPY pyproject.toml ./
COPY src ./src
COPY scripts/pc_sweep_qc.R ./pc_sweep_qc.R
RUN pip install --no-cache-dir .
ENTRYPOINT []
