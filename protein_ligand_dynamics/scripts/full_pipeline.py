#!/usr/bin/env python3
"""Полный пайплайн: от последовательности до анализа MD."""

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from analysis.trajectory_analyzer import TrajectoryAnalyzer
from docking.simple_docking import SimpleDocking
from md_simulation.mace_simulator import MACESimulator
from md_simulation.openmm_simulator import OpenMMSimulator
from structure_prediction.esmfold_predictor import ESMFoldPredictor
from utils.logging_config import setup_logger
from utils.visualization import MolecularVisualizer


class ProteinLigandPipeline:
    def __init__(self, config: dict, output_dir: str):
        self.config = config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = setup_logger(
            f"pipeline.{self.output_dir.resolve()}",
            self.output_dir / "pipeline.log",
        )
        self.results = {
            "start_time": datetime.now().isoformat(),
            "config": config,
            "steps": {},
        }

    def run(self) -> dict:
        self.logger.info("Запуск пайплайна")
        start_time = time.time()

        try:
            protein_pdb = self._step1_structure_prediction()
            complex_pdb = self._step2_docking(protein_pdb)
            trajectory_file = self._step3_md_simulation(complex_pdb)
            self._step4_analysis(trajectory_file, complex_pdb)
            self._step5_visualization(complex_pdb, trajectory_file)

            elapsed = (time.time() - start_time) / 60
            self.results["end_time"] = datetime.now().isoformat()
            self.results["elapsed_time_min"] = elapsed
            self.results["status"] = "success"
            self._save_results()
            self.logger.info("Пайплайн завершён за %.1f мин", elapsed)
            return self.results
        except Exception as exc:
            self.logger.exception("Ошибка пайплайна")
            self.results["status"] = "failed"
            self.results["error"] = str(exc)
            self._save_results()
            raise

    def _step1_structure_prediction(self) -> str:
        self.logger.info("Предсказание структуры белка")
        config = self.config["structure_prediction"]
        predictor = ESMFoldPredictor(use_local=config.get("use_local", False))
        output_pdb = self.output_dir / "protein.pdb"

        _, metadata = predictor.predict(
            sequence=config["sequence"],
            output_path=str(output_pdb),
        )
        self.results["steps"]["structure_prediction"] = {
            "sequence_length": len(config["sequence"]),
            "plddt_mean": metadata.get("plddt_mean", 0),
            "output_file": str(output_pdb),
        }
        return str(output_pdb)

    def _step2_docking(self, protein_pdb: str) -> str:
        self.logger.info("Докинг лиганда")
        config = self.config["docking"]
        docker = SimpleDocking(num_poses=config.get("num_poses", 10))
        docking_dir = self.output_dir / "docking"

        poses = docker.dock(
            protein_pdb=protein_pdb,
            ligand_smiles=config["ligand_smiles"],
            output_dir=str(docking_dir),
            binding_site_residues=config.get("binding_site"),
        )
        best_pose_pdb, best_energy = poses[0]

        complex_pdb = self.output_dir / "complex.pdb"
        docker.create_complex(protein_pdb, best_pose_pdb, str(complex_pdb))

        self.results["steps"]["docking"] = {
            "num_poses": len(poses),
            "best_energy": best_energy,
            "best_pose": best_pose_pdb,
            "complex_file": str(complex_pdb),
        }
        return str(complex_pdb)

    def _step3_md_simulation(self, complex_pdb: str) -> str:
        self.logger.info("Молекулярная динамика")
        config = self.config["md_simulation"]
        if config.get("use_mace", False):
            return self._run_mace_md(complex_pdb, config)
        return self._run_openmm_md(complex_pdb, config)

    def _run_mace_md(self, complex_pdb: str, config: dict) -> str:
        self.logger.info("Симуляция на MACE")
        simulator = MACESimulator(device=config.get("device", "cuda"))
        atoms = simulator.prepare_system(complex_pdb)
        atoms = simulator.minimize_energy(atoms, fmax=0.05)

        trajectory_file = self.output_dir / "trajectory_mace.dcd"
        md_results = simulator.run_md(
            atoms=atoms,
            temperature=config.get("temperature", 300.0),
            timestep=config.get("timestep", 1.0),
            num_steps=config.get("num_steps", 10000),
            trajectory_file=str(trajectory_file),
            log_interval=config.get("log_interval", 100),
        )
        self.results["steps"]["md_simulation"] = {"method": "MACE", **md_results}
        return str(trajectory_file)

    def _run_openmm_md(self, complex_pdb: str, config: dict) -> str:
        self.logger.info("Симуляция на OpenMM")
        simulator = OpenMMSimulator(
            complex_pdb=complex_pdb,
            use_gpu=config.get("use_gpu", True),
        )
        simulator.prepare(
            box_padding=config.get("box_padding", 1.0),
            ionic_strength=config.get("ionic_strength", 0.15),
        )
        simulator.equilibrate(
            npt_steps=config.get("equilibration_steps", 25000),
            nvt_steps=config.get("equilibration_steps", 25000),
            temperature=config.get("temperature", 300.0),
        )

        trajectory_file = self.output_dir / "trajectory.dcd"
        md_results = simulator.run_production(
            num_steps=config.get("num_steps", 500000),
            report_interval=config.get("report_interval", 1000),
            trajectory_file=str(trajectory_file),
        )
        self.results["steps"]["md_simulation"] = {"method": "OpenMM", **md_results}
        return str(trajectory_file)

    def _step4_analysis(self, trajectory_file: str, topology_file: str):
        self.logger.info("Анализ траектории")
        analyzer = TrajectoryAnalyzer(trajectory_file, topology_file)
        analysis_dir = self.output_dir / "analysis"
        results = analyzer.generate_report(str(analysis_dir))

        self.results["steps"]["analysis"] = {
            "rmsd_mean": float(np.mean(results["rmsd"])),
            "rmsd_std": float(np.std(results["rmsd"])),
            "rg_mean": float(np.mean(results["radius_of_gyration"])),
            "output_dir": str(analysis_dir),
        }

    def _step5_visualization(self, complex_pdb: str, trajectory_file: str):
        self.logger.info("Визуализация")
        MolecularVisualizer().plot_contact_map(
            complex_pdb,
            output_file=str(self.output_dir / "contact_map.png"),
        )

    def _save_results(self):
        results_file = self.output_dir / "results.json"
        with open(results_file, "w", encoding="utf-8") as handle:
            json.dump(
                self._json_safe(self.results),
                handle,
                indent=2,
                ensure_ascii=False,
            )
        self.logger.info("Результаты сохранены: %s", results_file)

    @staticmethod
    def _json_safe(obj):
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


def main():
    parser = argparse.ArgumentParser(description="Полный пайплайн: белок + лиганд + MD")
    parser.add_argument("--config", required=True, help="Путь к YAML-конфигу")
    parser.add_argument("--output", default="results", help="Каталог результатов")
    args = parser.parse_args()

    with open(args.config, encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    ProteinLigandPipeline(config, args.output).run()
    print(f"Результаты: {Path(args.output) / 'results.json'}")


if __name__ == "__main__":
    main()
