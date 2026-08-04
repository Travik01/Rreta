#!/usr/bin/env python3
"""
Полный пайплайн: от последовательности до анализа MD.

Использование:
    python scripts/01_full_pipeline.py --config config/pipeline_config.yaml
    python -m scripts.full_pipeline   --config config/pipeline_config.yaml
"""

import argparse
import yaml
import json
import time
from pathlib import Path
from datetime import datetime
import sys

import numpy as np

# Добавляем src в PATH
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from structure_prediction.esmfold_predictor import ESMFoldPredictor
from docking.simple_docking import SimpleDocking
from md_simulation.mace_simulator import MACESimulator
from md_simulation.openmm_simulator import OpenMMSimulator
from analysis.trajectory_analyzer import TrajectoryAnalyzer
from utils.visualization import MolecularVisualizer
from utils.logging_config import setup_logger


class ProteinLigandPipeline:
    """
    Полный production-ready пайплайн.

    Этапы:
    1. Предсказание структуры белка (ESMFold)
    2. Докинг лиганда (SimpleDocking/DiffDock)
    3. MD-симуляция (MACE или OpenMM)
    4. Анализ траектории
    5. Визуализация
    """

    def __init__(self, config: dict, output_dir: str):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Логирование
        self.logger = setup_logger(
            "Pipeline",
            self.output_dir / "pipeline.log"
        )

        # Результаты
        self.results = {
            "start_time": datetime.now().isoformat(),
            "config": config,
            "steps": {}
        }

    def run(self) -> dict:
        """Запуск полного пайплайна."""
        self.logger.info("="*70)
        self.logger.info("  ЗАПУСК ПАЙПЛАЙНА: Белок + Лиганд → MD → Анализ")
        self.logger.info("="*70)

        start_time = time.time()

        try:
            # Шаг 1: Структура белка
            protein_pdb = self._step1_structure_prediction()

            # Шаг 2: Докинг
            complex_pdb, best_pose_energy = self._step2_docking(protein_pdb)

            # Шаг 3: MD симуляция
            trajectory_file = self._step3_md_simulation(complex_pdb)

            # Шаг 4: Анализ
            analysis_results = self._step4_analysis(trajectory_file, complex_pdb)

            # Шаг 5: Визуализация
            self._step5_visualization(complex_pdb, trajectory_file)

            # Финализация
            self.results["end_time"] = datetime.now().isoformat()
            self.results["elapsed_time_min"] = (time.time() - start_time) / 60
            self.results["status"] = "success"

            self._save_results()

            self.logger.info("="*70)
            self.logger.info(f"  ✓ ПАЙПЛАЙН ЗАВЕРШЁН за {self.results['elapsed_time_min']:.1f} мин")
            self.logger.info(f"  Результаты: {self.output_dir}/")
            self.logger.info("="*70)

            return self.results

        except Exception as e:
            self.logger.error(f"ОШИБКА ПАЙПЛАЙНА: {e}", exc_info=True)
            self.results["status"] = "failed"
            self.results["error"] = str(e)
            self._save_results()
            raise

    def _step1_structure_prediction(self) -> str:
        """Шаг 1: Предсказание структуры белка."""
        self.logger.info("\n[ШАГ 1/5] Предсказание структуры белка")

        config = self.config["structure_prediction"]
        sequence = config["sequence"]

        predictor = ESMFoldPredictor(
            use_local=config.get("use_local", False)
        )

        output_pdb = self.output_dir / "protein.pdb"

        pdb_string, metadata = predictor.predict(
            sequence=sequence,
            output_path=str(output_pdb)
        )

        self.results["steps"]["structure_prediction"] = {
            "sequence_length": len(sequence),
            "plddt_mean": metadata.get("plddt_mean", 0),
            "output_file": str(output_pdb)
        }

        return str(output_pdb)

    def _step2_docking(self, protein_pdb: str) -> tuple:
        """Шаг 2: Докинг лиганда."""
        self.logger.info("\n[ШАГ 2/5] Молекулярный докинг")

        config = self.config["docking"]
        ligand_smiles = config["ligand_smiles"]

        docker = SimpleDocking(
            num_poses=config.get("num_poses", 10)
        )

        docking_dir = self.output_dir / "docking"

        poses = docker.dock(
            protein_pdb=protein_pdb,
            ligand_smiles=ligand_smiles,
            output_dir=str(docking_dir),
            binding_site_residues=config.get("binding_site", None)
        )

        # Выбор лучшей позы
        best_pose_pdb, best_energy = poses[0]

        # Создание комплекса
        complex_pdb = self.output_dir / "complex.pdb"
        docker.create_complex(
            protein_pdb,
            best_pose_pdb,
            str(complex_pdb)
        )

        self.results["steps"]["docking"] = {
            "num_poses": len(poses),
            "best_energy": best_energy,
            "best_pose": best_pose_pdb,
            "complex_file": str(complex_pdb)
        }

        return str(complex_pdb), best_energy

    def _step3_md_simulation(self, complex_pdb: str) -> str:
        """Шаг 3: MD-симуляция."""
        self.logger.info("\n[ШАГ 3/5] Молекулярная динамика")

        config = self.config["md_simulation"]
        use_mace = config.get("use_mace", True)

        if use_mace:
            return self._run_mace_md(complex_pdb, config)
        else:
            return self._run_openmm_md(complex_pdb, config)

    def _run_mace_md(self, complex_pdb: str, config: dict) -> str:
        """MD с MACE потенциалом."""
        self.logger.info("        Используем MACE ML-потенциал")

        simulator = MACESimulator(device=config.get("device", "cuda"))

        # Подготовка системы
        atoms = simulator.prepare_system(complex_pdb)

        # Минимизация
        atoms = simulator.minimize_energy(atoms, fmax=0.05)

        # Симуляция (DCD - читается анализатором)
        trajectory_file = self.output_dir / "trajectory_mace.dcd"

        md_results = simulator.run_md(
            atoms=atoms,
            temperature=config.get("temperature", 300.0),
            timestep=config.get("timestep", 1.0),
            num_steps=config.get("num_steps", 10000),
            trajectory_file=str(trajectory_file),
            log_interval=config.get("log_interval", 100)
        )

        self.results["steps"]["md_simulation"] = {
            "method": "MACE",
            **md_results
        }

        return str(trajectory_file)

    def _run_openmm_md(self, complex_pdb: str, config: dict) -> str:
        """MD с OpenMM (классическая)."""
        self.logger.info("        Используем OpenMM (AMBER силовое поле)")

        from md_simulation.openmm_simulator import OpenMMSimulator

        simulator = OpenMMSimulator(
            complex_pdb=complex_pdb,
            use_gpu=config.get("use_gpu", True)
        )

        # Подготовка
        simulator.prepare(
            box_padding=config.get("box_padding", 1.0),
            ionic_strength=config.get("ionic_strength", 0.15)
        )

        # Уравновешивание
        simulator.equilibrate(
            npt_steps=config.get("equilibration_steps", 25000),
            nvt_steps=config.get("equilibration_steps", 25000),
            temperature=config.get("temperature", 300.0)
        )

        # Production
        trajectory_file = self.output_dir / "trajectory.dcd"

        md_results = simulator.run_production(
            num_steps=config.get("num_steps", 500000),
            report_interval=config.get("report_interval", 1000),
            trajectory_file=str(trajectory_file)
        )

        self.results["steps"]["md_simulation"] = {
            "method": "OpenMM",
            **md_results
        }

        return str(trajectory_file)

    def _step4_analysis(self, trajectory_file: str, topology_file: str) -> dict:
        """Шаг 4: Анализ траектории."""
        self.logger.info("\n[ШАГ 4/5] Анализ траектории")

        analyzer = TrajectoryAnalyzer(
            trajectory_file=trajectory_file,
            topology_file=topology_file
        )

        analysis_dir = self.output_dir / "analysis"
        analysis_results = analyzer.generate_report(str(analysis_dir))

        self.results["steps"]["analysis"] = {
            "rmsd_mean": float(np.mean(analysis_results["rmsd"])),
            "rmsd_std": float(np.std(analysis_results["rmsd"])),
            "rg_mean": float(np.mean(analysis_results["radius_of_gyration"])),
            "output_dir": str(analysis_dir)
        }

        return analysis_results

    def _step5_visualization(self, complex_pdb: str, trajectory_file: str):
        """Шаг 5: Визуализация."""
        self.logger.info("\n[ШАГ 5/5] Визуализация")

        viz = MolecularVisualizer()

        # Карта контактов
        viz.plot_contact_map(
            complex_pdb,
            output_file=str(self.output_dir / "contact_map.png")
        )

        self.logger.info("        ✓ Визуализация завершена")

    def _save_results(self):
        """Сохранение финальных результатов."""
        results_file = self.output_dir / "results.json"

        with open(results_file, "w") as f:
            json.dump(self._json_safe(self.results), f, indent=2, ensure_ascii=False)

        self.logger.info(f"\n✓ Результаты сохранены: {results_file}")

    @staticmethod
    def _json_safe(obj):
        """Преобразование numpy-типов в обычные Python-типы для JSON."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, dict):
            return {str(k): ProteinLigandPipeline._json_safe(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [ProteinLigandPipeline._json_safe(v) for v in obj]
        return obj


# ============ ГЛАВНАЯ ФУНКЦИЯ ============

def main():
    parser = argparse.ArgumentParser(
        description="Полный пайплайн: Белок + Лиганд → MD"
    )

    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Путь к YAML конфигу"
    )

    parser.add_argument(
        "--output",
        type=str,
        default="results",
        help="Папка для результатов"
    )

    args = parser.parse_args()

    # Загрузка конфига
    with open(args.config) as f:
        config = yaml.safe_load(f)

    # Запуск пайплайна
    pipeline = ProteinLigandPipeline(config, args.output)
    results = pipeline.run()

    print("\n" + "="*70)
    print("  ФИНАЛЬНЫЕ РЕЗУЛЬТАТЫ")
    print("="*70)
    print(json.dumps(results["steps"], indent=2, default=str))


if __name__ == "__main__":
    main()
