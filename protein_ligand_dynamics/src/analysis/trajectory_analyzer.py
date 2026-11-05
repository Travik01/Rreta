"""Анализ MD-траекторий."""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

logger = logging.getLogger(__name__)


class TrajectoryAnalyzer:
    def __init__(self, trajectory_file: str, topology_file: str):
        self.traj = md.load(trajectory_file, top=topology_file)
        self.topology = self.traj.topology
        self.protein_indices = self.topology.select("protein")
        self.ligand_indices = self.topology.select("resname LIG or resname UNK")
        if len(self.ligand_indices) == 0:
            self.ligand_indices = self.topology.select("chainid 1")

        logger.info(
            "Загружено кадров: %d, атомов: %d",
            self.traj.n_frames,
            self.traj.n_atoms,
        )
        logger.info(
            "Атомов белка: %d, лиганда: %d",
            len(self.protein_indices),
            len(self.ligand_indices),
        )

    def compute_rmsd(
        self,
        reference_frame: int = 0,
        selection: str = "protein and name CA",
        align: bool = True,
    ) -> np.ndarray:
        atom_indices = self.topology.select(selection)
        if len(atom_indices) == 0:
            raise ValueError(f"Селекция не нашла атомов: {selection}")

        if align:
            self.traj.superpose(self.traj, frame=reference_frame, atom_indices=atom_indices)

        rmsd = md.rmsd(
            self.traj,
            self.traj,
            frame=reference_frame,
            atom_indices=atom_indices,
        ) * 10
        logger.info("RMSD: %.2f ± %.2f Å", np.mean(rmsd), np.std(rmsd))
        return rmsd

    def compute_rmsf(
        self,
        selection: str = "protein and name CA",
    ) -> tuple[np.ndarray, np.ndarray]:
        atom_indices = self.topology.select(selection)
        xyz = self.traj.xyz[:, atom_indices]
        mean_pos = xyz.mean(axis=0, keepdims=True)
        rmsf = np.sqrt(np.mean(np.sum((xyz - mean_pos) ** 2, axis=2), axis=0)) * 10
        residues = np.array([self.topology.atom(i).residue.resSeq for i in atom_indices])
        logger.info("RMSF: среднее %.2f Å", np.mean(rmsf))
        return rmsf, residues

    def compute_radius_of_gyration(self, selection: str = "protein") -> np.ndarray:
        atom_indices = self.topology.select(selection)
        rg = md.compute_rg(self.traj.atom_slice(atom_indices)) * 10
        logger.info("Rg: %.2f ± %.2f Å", np.mean(rg), np.std(rg))
        return rg

    def compute_contacts(self, cutoff: float = 0.45) -> tuple[np.ndarray, list]:
        if len(self.ligand_indices) == 0:
            logger.info("Лиганд не найден, контакты не рассчитаны")
            return np.array([]), []

        protein_atoms = np.asarray(self.protein_indices)
        column_of_atom = {index: col for col, index in enumerate(protein_atoms)}
        residue_of_atom = {}
        residue_ids = []
        for atom_index in protein_atoms:
            residue_id = self.topology.atom(atom_index).residue.index
            residue_of_atom[atom_index] = residue_id
            if residue_id not in residue_ids:
                residue_ids.append(residue_id)

        min_dist = np.full((len(self.traj), len(protein_atoms)), np.inf)
        for frame_index, frame in enumerate(self.traj):
            distances = cdist(
                frame.xyz[0, protein_atoms],
                frame.xyz[0, self.ligand_indices],
            )
            min_dist[frame_index] = distances.min(axis=1)

        contact_pairs = []
        frequencies = []
        for residue_id in residue_ids:
            columns = [
                column_of_atom[index]
                for index, rid in residue_of_atom.items()
                if rid == residue_id
            ]
            frequency = float(np.mean(np.any(min_dist[:, columns] < cutoff, axis=1)))
            if frequency > 0:
                contact_pairs.append((residue_id, "LIG"))
                frequencies.append(frequency)

        frequencies = np.asarray(frequencies)
        logger.info("Остатков в контакте с лигандом: %d", len(frequencies))
        return frequencies, contact_pairs

    def compute_ligand_protein_distance(self) -> np.ndarray:
        if len(self.ligand_indices) == 0:
            raise ValueError("Лиганд не найден")

        distances = []
        for frame in self.traj:
            matrix = cdist(
                frame.xyz[0, self.protein_indices],
                frame.xyz[0, self.ligand_indices],
            )
            distances.append(matrix.min() * 10)

        distances = np.asarray(distances)
        logger.info(
            "Расстояние белок–лиганд: %.2f ± %.2f Å",
            np.mean(distances),
            np.std(distances),
        )
        return distances

    def compute_hydrogen_bonds(
        self,
        distance_cutoff: float = 0.25,
        angle_cutoff: float = 120.0,
    ) -> list[tuple[int, int]]:
        if len(self.ligand_indices) == 0:
            return []

        hbonds = md.baker_hubbard(
            self.traj,
            distance_cutoff=distance_cutoff,
            angle_cutoff=angle_cutoff,
            periodic=False,
        )
        protein_set = set(self.protein_indices)
        ligand_set = set(self.ligand_indices)
        pairs = {
            (donor, acceptor)
            for donor, _, acceptor in hbonds
            if (donor in protein_set and acceptor in ligand_set)
            or (donor in ligand_set and acceptor in protein_set)
        }
        logger.info("Межмолекулярных H-связей: %d", len(pairs))
        return sorted(pairs)

    def generate_report(self, output_dir: str = "analysis_results") -> dict:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        results = {}

        rmsd = self.compute_rmsd()
        results["rmsd"] = rmsd
        self._plot_timeseries(rmsd, "RMSD (Å)", "RMSD", output_path / "rmsd.png")

        rmsf, residues = self.compute_rmsf()
        results["rmsf"] = rmsf
        results["rmsf_residues"] = residues
        self._plot_per_residue(
            residues,
            rmsf,
            "RMSF (Å)",
            "RMSF",
            output_path / "rmsf.png",
        )

        rg = self.compute_radius_of_gyration()
        results["radius_of_gyration"] = rg
        self._plot_timeseries(rg, "Rg (Å)", "Радиус гирации", output_path / "rg.png")

        if len(self.ligand_indices) > 0:
            contact_freq, contact_pairs = self.compute_contacts()
            results["contact_frequencies"] = contact_freq
            results["contact_pairs"] = contact_pairs

            distances = self.compute_ligand_protein_distance()
            results["ligand_distance"] = distances
            self._plot_timeseries(
                distances,
                "Расстояние (Å)",
                "Белок–лиганд",
                output_path / "distance.png",
            )

            results["hydrogen_bonds"] = self.compute_hydrogen_bonds()

        self._save_csv(results, output_path / "metrics.csv")
        logger.info("Отчёт сохранён: %s", output_path)
        return results

    @staticmethod
    def _plot_timeseries(data, ylabel, title, filename):
        plt.figure(figsize=(10, 4))
        plt.plot(data, linewidth=0.8)
        plt.xlabel("Кадр")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _plot_per_residue(residues, data, ylabel, title, filename):
        plt.figure(figsize=(12, 4))
        plt.bar(residues, data, width=1.0, alpha=0.7)
        plt.xlabel("Номер остатка")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches="tight")
        plt.close()

    @staticmethod
    def _save_csv(results, filename):
        data = {
            key: np.asarray(value, dtype=float)
            for key, value in results.items()
            if isinstance(value, np.ndarray) and value.ndim == 1 and len(value) > 0
        }

        if not data:
            return

        max_len = max(len(value) for value in data.values())
        aligned = {
            key: np.pad(value, (0, max_len - len(value)), constant_values=np.nan)
            for key, value in data.items()
        }
        pd.DataFrame(aligned).to_csv(filename, index=False)
        logger.info("Метрики сохранены: %s", filename)
