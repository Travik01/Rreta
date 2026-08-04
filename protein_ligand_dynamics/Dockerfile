FROM nvidia/cuda:11.8.0-cudnn8-runtime-ubuntu22.04

# Установка Miniconda
RUN apt-get update && apt-get install -y wget && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH=/opt/conda/bin:$PATH

# Копирование файлов
WORKDIR /app
COPY environment.yml .
COPY requirements.txt .

# Установка зависимостей
RUN conda env create -f environment.yml && \
    conda clean -afy

# Активация окружения
SHELL ["conda", "run", "-n", "protein_dynamics", "/bin/bash", "-c"]

# Копирование кода
COPY src/ ./src/
COPY scripts/ ./scripts/
COPY config/ ./config/

# Точка входа
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "protein_dynamics", "python", "scripts/01_full_pipeline.py"]
