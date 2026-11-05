FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends wget && \
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/conda/bin:$PATH

WORKDIR /app
COPY environment.yml .
RUN conda env create -f environment.yml && conda clean -afy

COPY src/ ./src/
COPY scripts/ ./scripts/
COPY config/ ./config/

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "protein_dynamics", "python", "scripts/01_full_pipeline.py"]
