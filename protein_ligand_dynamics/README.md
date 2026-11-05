# Protein-Ligand Molecular Dynamics Pipeline

A production-ready pipeline for predicting protein–ligand interactions through molecular dynamics simulation. The system integrates structure prediction, molecular docking, and trajectory analysis into a single automated workflow.

---

## Overview

The pipeline performs the following steps in sequence:

1. **3D structure prediction** of the target protein from its amino acid sequence (ESMFold or AlphaFold2).
2. **Molecular docking** of a small-molecule ligand into the predicted binding site (DiffDock, AutoDock Vina, or geometric docking).
3. **Molecular dynamics simulation** of the protein–ligand complex using either a machine-learning potential (MACE) or a classical force field (OpenMM).
4. **Trajectory analysis** — RMSD, RMSF, contact maps, hydrogen bonds, and binding-energy estimation.
5. **Visualization** — static plots, 3D structures, and animated trajectories.

---

## System Requirements

### Minimum (CPU-only)

| Component | Requirement |
|-----------|------------|
| Python | 3.10 or higher |
| RAM | 8 GB |
| Disk | 10 GB |

### Recommended (GPU-accelerated)

| Component | Requirement |
|-----------|------------|
| GPU | NVIDIA, 8+ GB VRAM |
| CUDA | 11.8 or higher |
| RAM | 16 GB |
| Disk | 50 GB |

---

## Installation

```bash
# Clone the repository
git clone https://github.com/your-repo/protein-ligand-md.git
cd protein-ligand-md

# Create the conda environment
conda env create -f environment.yml
conda activate protein_dynamics

# Verify installation
python -c "import mace; print('MACE OK')"
python -c "import openmm; print('OpenMM OK')"
```

---

## Quick Start

### Step 1 — Run the full pipeline

```bash
python scripts/01_full_pipeline.py \
    --config config/pipeline_config.yaml \
    --output my_first_run
```

### Step 2 — Analyze the trajectory

```bash
python scripts/02_analyze_results.py \
    --trajectory my_first_run/trajectory.dcd \
    --topology my_first_run/complex.pdb \
    --output my_first_run/analysis
```

### Step 3 — Visualize the results

```bash
python scripts/03_visualize.py \
    --pdb my_first_run/complex.pdb \
    --trajectory my_first_run/trajectory.dcd \
    --output my_first_run/viz \
    --create-gif
```

---

## Configuration

Edit `config/pipeline_config.yaml`:

```yaml
structure_prediction:
  sequence: "YOUR_PROTEIN_SEQUENCE"
  use_local: false        # set true for local ESMFold

docking:
  ligand_smiles: "YOUR_LIGAND_SMILES"
  num_poses: 20

md_simulation:
  use_mace: false         # set true for ML potential
  temperature: 300.0       # Kelvin
  num_steps: 500000        # ~1 ns at dt = 2 fs
```

---

## API Keys

All external APIs used in this pipeline are free of charge.

| Service | URL | API Key Required |
|---------|-----|:-----------------:|
| ESMFold | `https://api.esmatlas.com/foldSequence/v1/pdb/` | No |
| AlphaFold2 (Colab) | ColabFold notebook | Google account only |
| ChEMBL | `https://www.ebi.ac.uk/chembl/api/data/` | No |
| PubChem | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/` | No |

---

## Project Structure

```
protein_ligand_dynamics/
├── src/                         # Source code
│   ├── structure_prediction/     # ESMFold integration
│   ├── docking/                  # Docking modules
│   ├── md_simulation/            # MACE + OpenMM simulators
│   ├── analysis/                 # Trajectory analysis
│   └── utils/                    # Visualization, logging
│
├── scripts/                      # Executable scripts
│   ├── 01_full_pipeline.py      # Main pipeline
│   ├── 02_analyze_results.py     # Analysis
│   ├── 03_visualize.py           # Visualization
│   └── 04_batch_processing.py   # Batch processing
│
├── config/                      # Configuration files
├── tests/                        # Unit tests
├── notebooks/                    # Jupyter tutorials
├── data/examples/                # Example datasets
└── docs/                         # Documentation
```

---

## Usage Examples

### Example 1 — CDK2 kinase with a small-molecule inhibitor

```python
from scripts.full_pipeline import ProteinLigandPipeline

config = {
    "structure_prediction": {
        "sequence": "MENFQKVEKIGEGTYGVVYK...",
        "use_local": False
    },
    "docking": {
        "ligand_smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "num_poses": 10
    },
    "md_simulation": {
        "use_mace": False,
        "num_steps": 500000,
        "temperature": 310.0
    }
}

pipeline = ProteinLigandPipeline(config, "cdk2_results")
results = pipeline.run()
```

### Example 2 — Batch screening

Prepare a CSV file:

```csv
id,sequence,ligand_smiles
target1,MKFLK...,CC(=O)O
target2,AEQIV...,CN1C=NC2
target3,GSSGS...,C1CCCCC1
```

Run the batch processor:

```bash
python scripts/04_batch_processing.py \
    --input input.csv \
    --config config/pipeline_config.yaml \
    --output batch_results \
    --workers 4
```

---

## Result Interpretation

### RMSD (Root Mean Square Deviation)

| RMSD Range | Interpretation |
|------------|----------------|
| < 2 Å | Stable structure |
| 2–3 Å | Moderate conformational change |
| > 3 Å | Significant rearrangement |

### pLDDT (structure confidence score)

| pLDDT | Confidence |
|-------|------------|
| > 90 | Very high |
| 70–90 | Good |
| 50–70 | Low |
| < 50 | Unreliable — do not use |

### Binding energy (from docking)

| Energy (kcal/mol) | Binding strength |
|--------------------|-----------------|
| < −7 | Strong |
| −7 to −5 | Moderate |
| > −5 | Weak |

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| CUDA out of memory | Set `use_mace: false` to use OpenMM; reduce `num_steps` |
| ESMFold API timeout | Set `use_local: true` for a local ESMFold instance |
| OpenMM platform error | Verify CUDA with `nvidia-smi`; reinstall PyTorch with `pytorch-cuda=11.8` |

---

## References

- **ESMFold** — Lin et al., *Science*, 2023 — [GitHub](https://github.com/facebookresearch/esm)
- **MACE** — Batatia et al., *NeurIPS*, 2022 — [GitHub](https://github.com/ACEsuit/mace)
- **OpenMM** — Eastman et al., *JCTC*, 2017 — [Documentation](http://docs.openmm.org/)
- **MDAnalysis** — [Tutorials](https://www.mdanalysis.org/MDAnalysisTutorial/)

### Citation

```bibtex
@software{protein_ligand_md,
  title  = {Protein-Ligand MD Pipeline},
  author = {Your Name},
  year   = {2024},
  url    = {https://github.com/your-repo/protein-ligand-md}
}
```

---

## FAQ

**Is a GPU required?**
No. CPU mode is fully functional, but simulations run 10–100× slower.

**How long does a single run take?**
Structure prediction: ~30 min. Docking: ~10 min. MD simulation: 2–6 hours on GPU.

**Can I use custom force fields?**
Yes. Edit `src/md_simulation/openmm_simulator.py`.

**How do I add a new ligand?**
Provide its SMILES string in the config file. RDKit will generate the 3D conformation automatically.

---

## License

MIT License — free to use, modify, and distribute.

## Contributing

Pull requests are welcome. For substantial changes, please open an issue first to discuss the proposed modifications.

---
---

# Protein-Ligand Molecular Dynamics Pipeline

Production-готовый пайплайн для предсказания взаимодействия белков с лигандами методом молекулярной динамики. Система объединяет предсказание структуры, молекулярный докинг и анализ траекторий в единый автоматизированный workflow.

---

## Обзор

Пайплайн выполняет следующие этапы последовательно:

1. **Предсказание 3D-структуры** целевого белка по аминокислотной последовательности (ESMFold или AlphaFold2).
2. **Молекулярный докинг** лиганда в предсказанный сайт связывания (DiffDock, AutoDock Vina или геометрический докинг).
3. **Симуляция молекулярной динамики** комплекса белок–лиганд с использованием ML-потенциала (MACE) или классического силового поля (OpenMM).
4. **Анализ траекторий** — RMSD, RMSF, карты контактов, водородные связи, оценка энергии связывания.
5. **Визуализация** — статические графики, 3D-структуры и анимированные траектории.

---

## Системные требования

### Минимальные (CPU-only)

| Компонент | Требование |
|-----------|------------|
| Python | 3.10 и выше |
| ОЗУ | 8 ГБ |
| Диск | 10 ГБ |

### Рекомендуемые (GPU-ускорение)

| Компонент | Требование |
|-----------|------------|
| GPU | NVIDIA, от 8 ГБ VRAM |
| CUDA | 11.8 и выше |
| ОЗУ | 16 ГБ |
| Диск | 50 ГБ |

---

## Установка

```bash
# Клонирование репозитория
git clone https://github.com/your-repo/protein-ligand-md.git
cd protein-ligand-md

# Создание окружения
conda env create -f environment.yml
conda activate protein_dynamics

# Проверка установки
python -c "import mace; print('MACE OK')"
python -c "import openmm; print('OpenMM OK')"
```

---

## Быстрый старт

### Шаг 1 — Запуск полного пайплайна

```bash
python scripts/01_full_pipeline.py \
    --config config/pipeline_config.yaml \
    --output my_first_run
```

### Шаг 2 — Анализ траектории

```bash
python scripts/02_analyze_results.py \
    --trajectory my_first_run/trajectory.dcd \
    --topology my_first_run/complex.pdb \
    --output my_first_run/analysis
```

### Шаг 3 — Визуализация результатов

```bash
python scripts/03_visualize.py \
    --pdb my_first_run/complex.pdb \
    --trajectory my_first_run/trajectory.dcd \
    --output my_first_run/viz \
    --create-gif
```

---

## Конфигурация

Отредактируйте `config/pipeline_config.yaml`:

```yaml
structure_prediction:
  sequence: "YOUR_PROTEIN_SEQUENCE"
  use_local: false        # true для локального ESMFold

docking:
  ligand_smiles: "YOUR_LIGAND_SMILES"
  num_poses: 20

md_simulation:
  use_mace: false         # true для ML-потенциала
  temperature: 300.0       # Кельвин
  num_steps: 500000        # ~1 нс при dt = 2 фс
```

---

## API-ключи

Все внешние API, используемые в пайплайне, бесплатны.

| Сервис | URL | Нужен ключ |
|---------|-----|:----------:|
| ESMFold | `https://api.esmatlas.com/foldSequence/v1/pdb/` | Нет |
| AlphaFold2 (Colab) | ColabFold notebook | Только аккаунт Google |
| ChEMBL | `https://www.ebi.ac.uk/chembl/api/data/` | Нет |
| PubChem | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/` | Нет |

---

## Структура проекта

```
protein_ligand_dynamics/
├── src/                         # Исходный код
│   ├── structure_prediction/     # Интеграция ESMFold
│   ├── docking/                  # Модули докинга
│   ├── md_simulation/            # Симуляторы MACE + OpenMM
│   ├── analysis/                 # Анализ траекторий
│   └── utils/                    # Визуализация, логирование
│
├── scripts/                      # Исполняемые скрипты
│   ├── 01_full_pipeline.py      # Главный пайплайн
│   ├── 02_analyze_results.py     # Анализ
│   ├── 03_visualize.py           # Визуализация
│   └── 04_batch_processing.py   # Пакетная обработка
│
├── config/                      # Файлы конфигурации
├── tests/                        # Модульные тесты
├── notebooks/                    # Jupyter-туториалы
├── data/examples/                # Примеры данных
└── docs/                         # Документация
```

---

## Примеры использования

### Пример 1 — Киназа CDK2 с ингибитором

```python
from scripts.full_pipeline import ProteinLigandPipeline

config = {
    "structure_prediction": {
        "sequence": "MENFQKVEKIGEGTYGVVYK...",
        "use_local": False
    },
    "docking": {
        "ligand_smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "num_poses": 10
    },
    "md_simulation": {
        "use_mace": False,
        "num_steps": 500000,
        "temperature": 310.0
    }
}

pipeline = ProteinLigandPipeline(config, "cdk2_results")
results = pipeline.run()
```

### Пример 2 — Пакетный скрининг

Подготовьте CSV-файл:

```csv
id,sequence,ligand_smiles
target1,MKFLK...,CC(=O)O
target2,AEQIV...,CN1C=NC2
target3,GSSGS...,C1CCCCC1
```

Запуск пакетной обработки:

```bash
python scripts/04_batch_processing.py \
    --input input.csv \
    --config config/pipeline_config.yaml \
    --output batch_results \
    --workers 4
```

---

## Интерпретация результатов

### RMSD (среднеквадратичное отклонение)

| RMSD | Интерпретация |
|------|---------------|
| < 2 Å | Стабильная структура |
| 2–3 Å | Умеренные конформационные изменения |
| > 3 Å | Значительная перестройка |

### pLDDT (оценка достоверности структуры)

| pLDDT | Доверие |
|-------|---------|
| > 90 | Очень высокое |
| 70–90 | Хорошее |
| 50–70 | Низкое |
| < 50 | Недостоверно — не использовать |

### Энергия связывания (из докинга)

| Энергия (ккал/моль) | Сила связывания |
|---------------------|-----------------|
| < −7 | Сильное |
| −7 до −5 | Умеренное |
| > −5 | Слабое |

---

## Устранение неполадок

| Проблема | Решение |
|----------|---------|
| CUDA out of memory | Установите `use_mace: false` для OpenMM; уменьшите `num_steps` |
| ESMFold API timeout | Установите `use_local: true` для локального ESMFold |
| Ошибка платформы OpenMM | Проверьте CUDA командой `nvidia-smi`; переустановите PyTorch с `pytorch-cuda=11.8` |

---

## Ссылки

- **ESMFold** — Lin et al., *Science*, 2023 — [GitHub](https://github.com/facebookresearch/esm)
- **MACE** — Batatia et al., *NeurIPS*, 2022 — [GitHub](https://github.com/ACEsuit/mace)
- **OpenMM** — Eastman et al., *JCTC*, 2017 — [Документация](http://docs.openmm.org/)
- **MDAnalysis** — [Туториалы](https://www.mdanalysis.org/MDAnalysisTutorial/)

### Цитирование

```bibtex
@software{protein_ligand_md,
  title  = {Protein-Ligand MD Pipeline},
  author = {Maria Shaposhnik},
  year   = {2027},
  url    = {https://github.com/Travik01/Rreta/protein-ligand-md}
}
```

---

## FAQ

**Требуется ли GPU?**
Нет. CPU-режим полностью функционален, но симуляция работает в 10–100 раз медленнее.

**Сколько времени занимает один запуск?**
Предсказание структуры: ~30 мин. Докинг: ~10 мин. MD-симуляция: 2–6 часов на GPU.

**Можно ли использовать свои силовые поля?**
Да. Отредактируйте `src/md_simulation/openmm_simulator.py`.

**Как добавить новый лиганд?**
Укажите его SMILES-строку в конфигурации. RDKit автоматически сгенерирует 3D-конформацию.

---

## Лицензия

MIT License — свободное использование, модификация и распространение.

## Контрибьюты

Pull requests приветствуются. Для существенных изменений откройте issue для предварительного обсуждения.
