# Установка

## Conda

```bash
conda env create -f environment.yml
conda activate protein_dynamics
python -c "import mdtraj, rdkit, openmm; print('ok')"
```

Для CPU без CUDA удалите строку `pytorch::pytorch-cuda=11.8` из `environment.yml` перед созданием окружения.

## pip (CPU)

```bash
python3.10 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

`requirements.txt` не содержит PyTorch, MACE и fair-esm. Их устанавливают отдельно, если нужны MACE или локальный ESMFold:

```bash
pip install torch mace-torch fair-esm
```

RDKit и OpenMM надёжнее ставить из conda-forge; при ошибках установки через pip используйте вариант с Conda.

## Docker

```bash
docker build -t protein-ligand-md .
docker run -v "$(pwd)/results:/app/results" protein-ligand-md \
    --config config/pipeline_config.yaml --output results
```

## Проверка

```bash
python --version
python -c "import torch; print(torch.cuda.is_available())"
python -c "import mdtraj, rdkit, openmm; print('ok')"
pytest
```

## Дополнительные компоненты

- `py3Dmol` — интерактивный просмотр структур в Jupyter (`MolecularVisualizer.show_structure`).
- `fair-esm` — локальный ESMFold, нужен PyTorch и большой объём видеопамяти.
- `mace-torch` — симуляции с потенциалом MACE, нужен PyTorch.

## Частые проблемы

| Проблема | Что проверить |
|----------|---------------|
| CUDA не видна | `nvidia-smi`, версию драйвера и сборку PyTorch |
| OpenMM не импортируется | Установить пакет из conda-forge |
| Ошибки RDKit при pip | Использовать `conda install -c conda-forge rdkit` |
| ESMFold API не отвечает | Доступ в интернет, длину последовательности, повторить запрос |
