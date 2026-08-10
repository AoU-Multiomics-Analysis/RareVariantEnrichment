FROM python:3.12-slim-bookworm
RUN apt-get update \
    && apt-get install -y --no-install-recommends tabix ca-certificates \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/rare-variant-enrichment
COPY pyproject.toml ./
COPY src ./src
RUN pip install --no-cache-dir .
ENTRYPOINT []
