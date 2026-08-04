#!/usr/bin/env python3
"""
Пакетная обработка множества белков/лигандов.

Использование:
    python scripts/04_batch_processing.py --input batch_input.csv \
                                          --output batch_results \
                                          --config config/pipeline_config.yaml
"""

import argparse
import pandas as pd
import yaml
import json
from pathlib import Path
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))

from full_pipeline import ProteinLigandPipeline


def process_single(row: dict, base_config: dict, output_base: Path, job_id: str) -> dict:
    """Обработка одного белка/лиганда."""
    
    # Обновление конфига
    config = base_config.copy()
    config["structure_prediction"]["sequence"] = row["sequence"]
    config["docking"]["ligand_smiles"] = row["ligand_smiles"]
    
    # Выходная папка для этой задачи
    output_dir = output_base / job_id
    
    try:
        pipeline = ProteinLigandPipeline(config, str(output_dir))
        results = pipeline.run()
        
        return {
            "id": job_id,
            "status": "success",
            "rmsd_mean": results["steps"]["analysis"]["rmsd_mean"],
            "output_dir": str(output_dir)
        }
    
    except Exception as e:
        return {
            "id": job_id,
            "status": "failed",
            "error": str(e)
        }


def main():
    parser = argparse.ArgumentParser(description="Пакетная обработка")
    
    parser.add_argument(
        "--input",
        required=True,
        help="CSV файл с колонками: id, sequence, ligand_smiles"
    )
    
    parser.add_argument(
        "--output",
        default="batch_results",
        help="Базовая папка для результатов"
    )
    
    parser.add_argument(
        "--config",
        required=True,
        help="Базовый YAML конфиг"
    )
    
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Количество параллельных задач"
    )
    
    args = parser.parse_args()
    
    # Загрузка данных
    df = pd.read_csv(args.input)
    print(f"[Batch] Загружено {len(df)} задач")
    
    # Проверка колонок
    required_cols = ["sequence", "ligand_smiles"]
    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(f"Отсутствуют колонки: {missing}")
    
    # Загрузка базового конфига
    with open(args.config) as f:
        base_config = yaml.safe_load(f)
    
    output_base = Path(args.output)
    output_base.mkdir(parents=True, exist_ok=True)
    
    # Обработка
    results = []
    
    # Формирование job_id для каждой строки (берём колонку "id" или индекс)
    def _job_id(idx, row):
        if "id" in df.columns and pd.notna(row.get("id")):
            return str(row["id"])
        return f"job_{idx}"

    if args.workers == 1:
        # Последовательная обработка
        for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
            result = process_single(
                row.to_dict(), base_config, output_base, _job_id(idx, row)
            )
            results.append(result)
    
    else:
        # Параллельная обработка
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    process_single, row.to_dict(), base_config, output_base, _job_id(idx, row)
                ): idx
                for idx, row in df.iterrows()
            }
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Processing"):
                result = future.result()
                results.append(result)
    
    # Сохранение результатов
    results_df = pd.DataFrame(results)
    results_file = output_base / "batch_results.csv"
    results_df.to_csv(results_file, index=False)
    
    # Сводка
    success_count = (results_df["status"] == "success").sum()
    failed_count = (results_df["status"] == "failed").sum()
    
    print("\n" + "="*70)
    print("  ПАКЕТНАЯ ОБРАБОТКА ЗАВЕРШЕНА")
    print("="*70)
    print(f"  Успешно: {success_count}/{len(df)}")
    print(f"  Ошибки: {failed_count}/{len(df)}")
    print(f"  Результаты: {results_file}")
    print("="*70)


if __name__ == "__main__":
    main()
