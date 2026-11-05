# Protein-Ligand Dynamics

Пайплайн для подготовки структуры белка, докинга лиганда, молекулярной динамики и анализа траектории.

## Что делает

1. Получает PDB по аминокислотной последовательности через ESMFold (публичный API или локальная модель).
2. Строит 3D-конформер лиганда из SMILES и размещает его около центра белка или указанных остатков (RDKit).
3. Запускает MD на OpenMM (по умолчанию) или MACE.
4. Считает RMSD, RMSF, радиус гирации, контакты и минимальное расстояние белок–лиганд (MDTraj).
5. Сохраняет графики, `metrics.csv`, `summary.json` и `results.json`.

## Требования

- Python 3.10+
- Conda (основной вариант) либо pip для CPU-режима
- Для локального ESMFold и MACE: PyTorch с CUDA и соответствующие пакеты
- Для ESMFold API: доступ в интернет

## Установка

```bash
conda env create -f environment.yml
conda activate protein_dynamics
```

Установка через pip и запуск в Docker описаны в [docs/installation.md](docs/installation.md).

## Запуск

```bash
python scripts/01_full_pipeline.py \
    --config config/pipeline_config.yaml \
    --output results
```

Отдельные этапы:

```bash
python scripts/02_analyze_results.py \
    --trajectory results/trajectory.dcd \
    --topology results/complex.pdb \
    --output results/analysis

python scripts/03_visualize.py \
    --pdb results/complex.pdb \
    --trajectory results/trajectory.dcd \
    --output results/viz \
    --create-gif
```

Пакетная обработка CSV с колонками `id`, `sequence`, `ligand_smiles`:

```bash
python scripts/04_batch_processing.py \
    --input input.csv \
    --config config/pipeline_config.yaml \
    --output batch_results \
    --workers 4
```

## Конфигурация

`config/pipeline_config.yaml`:

| Ключ | Назначение |
|------|------------|
| `structure_prediction.sequence` | Аминокислотная последовательность |
| `structure_prediction.use_local` | `false` — ESMFold API, `true` — локальная модель |
| `docking.ligand_smiles` | SMILES лиганда |
| `docking.num_poses` | Число генерируемых поз |
| `docking.binding_site` | Список номеров остатков; `null` — центр белка |
| `md_simulation.use_mace` | `false` — OpenMM, `true` — MACE |
| `md_simulation.temperature` | Температура, K |
| `md_simulation.num_steps` | Число шагов MD |
| `md_simulation.device` | `cuda` или `cpu` для MACE |
| `md_simulation.use_gpu` | Платформа OpenMM |
| `md_simulation.timestep`, `log_interval` | Параметры MACE |
| `md_simulation.box_padding`, `ionic_strength`, `equilibration_steps`, `report_interval` | Параметры OpenMM |

## Результаты

| Файл | Содержимое |
|------|------------|
| `protein.pdb` | Структура белка |
| `complex.pdb` | Белок с лучшей позой лиганда |
| `trajectory.dcd` / `trajectory_mace.dcd` | Траектория |
| `analysis/` | Графики и `metrics.csv` |
| `contact_map.png` | Карта контактов C-alpha |
| `results.json` | Сводка по этапам |
| `pipeline.log` | Лог запуска |

## Тесты

```bash
pytest
```

Медленные тесты (сеть, GPU, полный пайплайн) отключены по умолчанию через `pytest.ini`.

## Структура

```text
protein_ligand_dynamics/
├── config/pipeline_config.yaml
├── docs/installation.md
├── scripts/               # CLI-скрипты и точки входа с номерами
├── src/
│   ├── analysis/
│   ├── docking/
│   ├── md_simulation/
│   ├── structure_prediction/
│   └── utils/
└── tests/
```

## Ограничения

- Докинг геометрический: он не оценивает аффинность связывания. Значение `energy` у позы — энергия силового поля лиганда, а не энергия связывания с белком.
- OpenMM использует `amber14-all.xml`. Для корректной симуляции комплекса с лигандом нужны параметры лиганда и проверка подготовленной системы под конкретную задачу.
- ESMFold API — внешний сервис; возможны задержки и ограничения на длину последовательности.
- MACE и локальный ESMFold требуют дополнительных зависимостей и GPU.

## Ссылки

- [ESMFold](https://github.com/facebookresearch/esm)
- [OpenMM](https://openmm.org/)
- [MACE](https://github.com/ACEsuit/mace)
- [MDTraj](https://www.mdtraj.org/)
- [RDKit](https://www.rdkit.org/)
