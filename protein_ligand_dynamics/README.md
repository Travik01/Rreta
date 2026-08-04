# Protein-Ligand Molecular Dynamics Pipeline

🔬 **Production-ready пайплайн для предсказания взаимодействия белков с лигандами через молекулярную динамику.**

## 🌟 Возможности

- ✅ **Предсказание 3D структуры белка** (ESMFold/AlphaFold2)
- ✅ **Молекулярный докинг** (DiffDock/AutoDock Vina/геометрический)
- ✅ **MD-симуляция** (MACE ML-потенциал или OpenMM классическая)
- ✅ **Комплексный анализ траекторий** (RMSD, RMSF, контакты, H-связи)
- ✅ **Визуализация** (3D структуры, графики, анимации)
- ✅ **Пакетная обработка** (multiple targets)

## 📋 Требования

### Минимальные (CPU-only)
- Python 3.10+
- 8 GB RAM
- 10 GB диск

### Рекомендуемые (GPU)
- NVIDIA GPU (8+ GB VRAM)
- CUDA 11.8+
- 16 GB RAM
- 50 GB диск

## 🚀 Быстрый старт

### 1. Установка

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

### 2. Базовый пример

```bash
# Запуск с примером конфига
python scripts/01_full_pipeline.py \
    --config config/pipeline_config.yaml \
    --output my_first_run
```

### 3. Анализ результатов

```bash
# Анализ траектории
python scripts/02_analyze_results.py \
    --trajectory my_first_run/trajectory.dcd \
    --topology my_first_run/complex.pdb \
    --output my_first_run/analysis

# Визуализация
python scripts/03_visualize.py \
    --pdb my_first_run/complex.pdb \
    --trajectory my_first_run/trajectory.dcd \
    --output my_first_run/viz \
    --create-gif
```

## 📁 Структура проекта

```
protein_ligand_dynamics/
├── src/                        # Исходный код
│   ├── structure_prediction/   # ESMFold предсказание
│   ├── docking/                # Докинг модули
│   ├── md_simulation/          # MACE + OpenMM симуляторы
│   ├── analysis/               # Анализ траекторий
│   └── utils/                  # Визуализация, логирование
│
├── scripts/                    # Исполняемые скрипты
│   ├── 01_full_pipeline.py     # Главный пайплайн
│   ├── 02_analyze_results.py   # Анализ
│   ├── 03_visualize.py         # Визуализация
│   └── 04_batch_processing.py  # Пакетная обработка
│
├── config/                     # Конфигурации
│   ├── pipeline_config.yaml    # Основной конфиг
│   └── api_keys.yaml.template  # Шаблон для API ключей
│
├── tests/                      # Тесты
├── notebooks/                  # Jupyter туториалы
├── data/examples/              # Примеры данных
└── docs/                       # Документация
```

## 🔧 Конфигурация

Отредактируйте `config/pipeline_config.yaml`:

```yaml
structure_prediction:
  sequence: "YOUR_PROTEIN_SEQUENCE"
  use_local: false  # true для локального ESMFold

docking:
  ligand_smiles: "YOUR_LIGAND_SMILES"
  num_poses: 20

md_simulation:
  use_mace: false  # true для ML-потенциала
  temperature: 300.0
  num_steps: 500000  # ~1 нс при dt=2fs
```

## 🔑 API Ключи (все БЕСПЛАТНЫЕ!)

### ESMFold API
```yaml
URL: https://api.esmatlas.com/foldSequence/v1/pdb/
Ключ: НЕ ТРЕБУЕТСЯ ✅
Лимит: До 800 а.о., без rate limit
```

### AlphaFold2 (через Google Colab)
```yaml
URL: https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb
Ключ: Google Account
GPU: Бесплатно 12 часов/день ✅
```

### ChEMBL (химические данные)
```yaml
URL: https://www.ebi.ac.uk/chembl/api/data/
Ключ: НЕ ТРЕБУЕТСЯ ✅
```

### PubChem
```yaml
URL: https://pubchem.ncbi.nlm.nih.gov/rest/pug/
Ключ: НЕ ТРЕБУЕТСЯ ✅
```

## 📊 Примеры использования

### Пример 1: CDK2 киназа + ингибитор

```python
from scripts.full_pipeline import ProteinLigandPipeline

config = {
    "structure_prediction": {
        "sequence": "MENFQKVEKIGEGTYGVVYK...",  # CDK2 sequence
        "use_local": False
    },
    "docking": {
        "ligand_smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "num_poses": 10
    },
    "md_simulation": {
        "use_mace": False,
        "num_steps": 500000,
        "temperature": 310.0  # 37°C
    }
}

pipeline = ProteinLigandPipeline(config, "cdk2_results")
results = pipeline.run()
```

### Пример 2: Пакетный скрининг

```csv
# input.csv
id,sequence,ligand_smiles
target1,MKFLK...,CC(=O)O
target2,AEQIV...,CN1C=NC2
target3,GSSGS...,C1CCCCC1
```

```bash
python scripts/04_batch_processing.py \
    --input input.csv \
    --config config/pipeline_config.yaml \
    --output batch_results \
    --workers 4
```

## 📈 Интерпретация результатов

### RMSD (Root Mean Square Deviation)
- **< 2 Å**: Стабильная структура ✅
- **2-3 Å**: Умеренные изменения ⚠️
- **> 3 Å**: Значительная переориентация ❌

### pLDDT (confidence score структуры)
- **> 90**: Очень высокая уверенность ✅
- **70-90**: Хорошая модель ✅
- **50-70**: Низкая уверенность ⚠️
- **< 50**: Не используй ❌

### Энергия связывания (из докинга)
- **< -7 ккал/моль**: Сильное связывание
- **-7 до -5**: Умеренное
- **> -5**: Слабое

## 🐛 Troubleshooting

### Проблема: "CUDA out of memory"
```bash
# Решение 1: Уменьши размер системы
use_mace: false  # Используй OpenMM

# Решение 2: Уменьши batch size
md_simulation:
  num_steps: 100000  # Вместо 500000
```

### Проблема: "ESMFold API timeout"
```bash
# Решение: Используй локальную версию
structure_prediction:
  use_local: true
```

### Проблема: "OpenMM platform error"
```bash
# Проверь CUDA
nvidia-smi

# Переустанови с правильной версией
conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
```

## 📚 Дополнительные ресурсы

- [ESMFold документация](https://github.com/facebookresearch/esm)
- [MACE GitHub](https://github.com/ACEsuit/mace)
- [OpenMM документация](http://docs.openmm.org/)
- [MDAnalysis туториалы](https://www.mdanalysis.org/MDAnalysisTutorial/)

## 📝 Цитирование

Если используешь этот пайплайн, пожалуйста, цитируй:

```bibtex
@software{protein_ligand_md,
  title={Protein-Ligand MD Pipeline},
  author={Your Name},
  year={2024},
  url={https://github.com/your-repo/protein-ligand-md}
}
```

И оригинальные методы:
- **ESMFold**: Lin et al., Science 2023
- **MACE**: Batatia et al., NeurIPS 2022
- **OpenMM**: Eastman et al., JCTC 2017

## 📄 Лицензия

MIT License - свободно используй в своих проектах!

## 🤝 Контрибьюции

Pull requests приветствуются! Для больших изменений открывайте issue.

## ❓ FAQ

**Q: Нужен ли мне GPU?**  
A: Нет, но рекомендуется. CPU режим работает, но в 10-100x медленнее.

**Q: Сколько времени занимает один запуск?**  
A: ~30 мин (структура) + ~10 мин (докинг) + ~2-6 часов (MD на GPU).

**Q: Можно ли использовать свои силовые поля?**  
A: Да! Отредактируй `src/md_simulation/openmm_simulator.py`.

**Q: Как добавить новый лиганд?**  
A: Просто укажи SMILES в конфиге. RDKit автоматически сгенерирует 3D.

---

**💡 Нашли баг? Есть вопрос?** Открывайте [issue](https://github.com/your-repo/issues)!
```

---
