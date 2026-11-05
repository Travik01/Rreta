#!/usr/bin/env python3
# ruff: noqa: N999
"""Пакетная обработка белков и лигандов."""

import argparse
import copy
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))

from full_pipeline import ProteinLigandPipeline


def process_single(row: dict, base_config: dict, output_base: Path, job_id: str) -> dict:
    config = copy.deepcopy(base_config)
    config["structure_prediction"]["sequence"] = row["sequence"]
    config["docking"]["ligand_smiles"] = row["ligand_smiles"]
    output_dir = output_base / job_id

    try:
        results = ProteinLigandPipeline(config, str(output_dir)).run()
        return {
            "id": job_id,
            "status": "success",
            "rmsd_mean": results["steps"]["analysis"]["rmsd_mean"],
            "output_dir": str(output_dir),
        }
    except Exception as exc:  # noqa: BLE001 - сбой задачи фиксируется в отчёте
        return {"id": job_id, "status": "failed", "error": str(exc)}


def main():
    parser = argparse.ArgumentParser(description="Пакетная обработка")
    parser.add_argument(
        "--input",
        required=True,
        help="CSV с колонками id, sequence, ligand_smiles",
    )
    parser.add_argument("--output", default="batch_results", help="Каталог результатов")
    parser.add_argument("--config", required=True, help="Базовый YAML-конфиг")
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Число параллельных задач",
    )
    args = parser.parse_args()

    if args.workers < 1:
        parser.error("--workers должен быть не меньше 1")

    df = pd.read_csv(args.input)
    missing = {"sequence", "ligand_smiles"} - set(df.columns)
    if missing:
        raise ValueError(f"Отсутствуют колонки: {sorted(missing)}")

    with open(args.config, encoding="utf-8") as handle:
        base_config = yaml.safe_load(handle)

    output_base = Path(args.output)
    output_base.mkdir(parents=True, exist_ok=True)

    def job_id(idx, row):
        if "id" in df.columns and pd.notna(row.get("id")):
            return str(row["id"])
        return f"job_{idx}"

    jobs = [(row.to_dict(), job_id(idx, row)) for idx, row in df.iterrows()]

    results = []
    if args.workers == 1:
        for row, jid in tqdm(jobs, desc="Processing"):
            results.append(process_single(row, base_config, output_base, jid))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [
                executor.submit(process_single, row, base_config, output_base, jid)
                for row, jid in jobs
            ]
            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Processing"
            ):
                results.append(future.result())

    results_df = pd.DataFrame(results)
    results_file = output_base / "batch_results.csv"
    results_df.to_csv(results_file, index=False)

    success = int((results_df["status"] == "success").sum())
    failed = int((results_df["status"] == "failed").sum())
    print(f"Успешно: {success}/{len(df)}")
    print(f"Ошибки: {failed}/{len(df)}")
    print(f"Результаты: {results_file}")


if __name__ == "__main__":
    main()
