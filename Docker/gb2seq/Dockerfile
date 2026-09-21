# syntax=docker/dockerfile:1

# ============================================================
# 1. Common Conda base
# ============================================================
FROM mambaorg/micromamba:1.5.10 AS base

RUN micromamba install -y -n base \
        -c conda-forge \
        -c bioconda \
        --strict-channel-priority \
        mafft=7.526 \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

SHELL ["/bin/bash", "-lc"]


# ============================================================
# 2. Runtime image used by Nextflow
# ============================================================
FROM base AS runtime

ARG GB2SEQ_COMMIT

LABEL org.opencontainers.image.source="https://github.com/<owner>/gb2seq"
LABEL org.opencontainers.image.revision="${GB2SEQ_COMMIT}"
LABEL org.opencontainers.image.description="gb2seq execution environment with MAFFT and pandas"

COPY --from=ghcr.io/astral-sh/uv:0.12.15 \
    /uv /uvx /usr/local/bin/

ENV UV_PROJECT_ENVIRONMENT=/opt/conda/envs/gb2seq-venv
ENV UV_PYTHON=3.10.21

ENV PATH=/opt/gb2seq/bin:/opt/conda/envs/gb2seq-venv/bin:/opt/conda/bin:$PATH

RUN mkdir -p /opt/conda/envs/gb2seq-venv && \
    chown -R $MAMBA_USER:$MAMBA_USER /opt/conda/envs/gb2seq-venv

# Everything copied here belongs to the last fixed gb2seq commit.
WORKDIR /opt/gb2seq

COPY pyproject.toml uv.lock README.md .python-version ./
COPY bin ./bin
COPY src ./src

# Install gb2seq and all locked runtime dependencies in one operation.
# --no-editable ensures that imports do not depend on /opt/gb2seq/src.
RUN uv sync \
    --frozen \
    --no-dev \
    --no-editable

# pandas is required by the nextflow pipeline scripts, although it is not a
# dependency of gb2seq itself.
RUN uv pip install \
    --python /opt/conda/envs/gb2seq-venv/bin/python \
    pandas==2.3.3

WORKDIR /workspace


# ============================================================
# 3. Development and unit-test image
# ============================================================
FROM runtime AS development

RUN uv pip install \
    --python /opt/conda/envs/gb2seq-venv/bin/python \
    pytest==9.1.1 \
    jupyterlab==4.6.3 \
    ipykernel==7.3.0

CMD [\
    "jupyter",\
    "lab",\
    "--ip=0.0.0.0",\
    "--no-browser"\
]