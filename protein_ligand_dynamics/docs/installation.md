# Подробная инструкция по установке

## Вариант 1: Conda (рекомендуется)

### Шаг 1: Установка Miniconda

**Linux/Mac:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

**Windows:**
Скачай с https://docs.conda.io/en/latest/miniconda.html

### Шаг 2: Создание окружения

```bash
conda env create -f environment.yml
conda activate protein_dynamics
```

### Шаг 3: Проверка

```bash
# Python версия
python --version  # Должно быть 3.10+

# CUDA (если есть GPU)
python -c "import torch; print(torch.cuda.is_available())"

# Ключевые пакеты
python -c "import mace; import openmm; import mdtraj; print('OK')"
```

## Вариант 2: pip (только CPU)

```bash
# Создание venv
python3.10 -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Установка
pip install -r requirements.txt
```

## Вариант 3: Docker

```bash
# Сборка
docker build -t protein-md .

# Запуск
docker run -v $(pwd)/data:/data protein-md --config /data/config.yaml
```

## Опциональные компоненты

### DiffDock (продвинутый докинг)

```bash
git clone https://github.com/gcorso/DiffDock.git external/DiffDock
cd external/DiffDock
pip install -e .
```

### ESMFold локально (требует GPU 16GB+)

```bash
pip install fair-esm
```

### AutoDock Vina

```bash
conda install -c conda-forge vina
```

## Troubleshooting установки

### CUDA не найдена
```bash
# Проверь версию драйвера
nvidia-smi

# Переустанови PyTorch с правильной CUDA
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
```

### OpenMM ошибка
```bash
# Переустанови
conda remove openmm
conda install -c conda-forge openmm
```

### RDKit ошибка
```bash
# Используй conda вместо pip
conda install -c conda-forge rdkit
```
```

---

## 🎯 Полный план запуска (пошаговая инструкция)

### Шаг 1: Подготовка окружения (15 мин)

```bash
# 1. Клонирование (если репозиторий создан)
git clone <URL>
cd protein_ligand_dynamics

# 2. Создание окружения
conda env create -f environment.yml

# 3. Активация
conda activate protein_dynamics

# 4. Проверка
python -c "import torch; print('PyTorch:', torch.__version__)"
python -c "import mace; print('MACE: OK')"
python -c "import openmm; print('OpenMM: OK')"
```

### Шаг 2: Тестовый запуск (5 мин)

```bash
# Быстрый тест на маленьком белке
python tests/test_structure_prediction.py
python tests/test_docking.py
```

### Шаг 3: Подготовка данных

Создай файл `my_project/input.yaml`:

```yaml
structure_prediction:
  sequence: |
    MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNH
    PNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHS
  use_local: false

docking:
  ligand_smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
  num_poses: 10

md_simulation:
  use_mace: false
  use_gpu: true
  temperature: 300.0
  num_steps: 100000
  report_interval: 1000
```

### Шаг 4: Запуск пайплайна (2-4 часа)

```bash
python scripts/01_full_pipeline.py \
    --config my_project/input.yaml \
    --output my_project/results
```

### Шаг 5: Анализ (5 мин)

```bash
python scripts/02_analyze_results.py \
    --trajectory my_project/results/trajectory.dcd \
    --topology my_project/results/complex.pdb \
    --output my_project/results/analysis
```

### Шаг 6: Визуализация

```bash
python scripts/03_visualize.py \
    --pdb my_project/results/complex.pdb \
    --trajectory my_project/results/trajectory.dcd \
    --output my_project/results/viz \
    --create-gif
```

### Шаг 7: Просмотр результатов

```bash
# Смотри метрики
cat my_project/results/analysis/summary.json

# Графики
open my_project/results/analysis/combined_analysis.png
open my_project/results/viz/trajectory.gif
```

---

**Готово!** 🎉 Полный production-ready пайплайн с MACE, OpenMM, всеми скриптами, тестами и документацией!
