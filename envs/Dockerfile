FROM mambaorg/micromamba:2.8.1

USER root
RUN mkdir -p /opt/rare-variant-enrichment \
    && chown -R "${MAMBA_USER}:${MAMBA_USER}" /opt/rare-variant-enrichment

USER $MAMBA_USER
RUN micromamba install --yes --name base --strict-channel-priority \
        -c conda-forge \
        -c bioconda \
        python=3.12.14 \
        numpy=2.5.2 \
        htslib=1.24 \
        r-base=4.6.1 \
        r-tidyverse=2.0.0 \
        r-ggplot2=4.0.3 \
        r-ggrepel=0.9.8 \
        pip=26.2.1 \
        setuptools=84.0.0 \
        wheel=0.48.0 \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/bin:${PATH}"
WORKDIR /opt/rare-variant-enrichment
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml ./
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src
COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/pc_sweep_qc.R ./pc_sweep_qc.R
RUN /opt/conda/bin/python -m pip install --no-cache-dir --no-deps --no-build-isolation .

ENTRYPOINT []
