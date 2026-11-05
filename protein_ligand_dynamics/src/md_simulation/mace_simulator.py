"""Молекулярная динамика с потенциалом MACE."""

import logging
from pathlib import Path

import numpy as np
from ase import Atoms, units
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

try:
    import torch
except ImportError:
    torch = None

try:
    from mace.calculators import MACECalculator
except ImportError:
    MACECalculator = None

logger = logging.getLogger(__name__)


class MACESimulator:
    def __init__(
        self,
        model_path: str = "mace_mp_0",
        device: str = "cuda",
        default_dtype: str = "float32",
    ):
        if torch is None or MACECalculator is None:
            raise ImportError(
                "Для MACE установите PyTorch и mace-torch: pip install torch mace-torch"
            )

        if device == "cuda" and not torch.cuda.is_available():
            device = "cpu"
        self.device = device

        logger.info("Инициализация MACE на %s", self.device)
        try:
            self.calculator = MACECalculator(
                model_paths=model_path,
                device=self.device,
                default_dtype=default_dtype,
            )
        except Exception as exc:
            raise RuntimeError(
                f"Не удалось загрузить модель MACE '{model_path}': {exc}"
            ) from exc

    def prepare_system(self, pdb_file: str, remove_hydrogens: bool = False) -> Atoms:
        atoms = ase_read(pdb_file)
        if remove_hydrogens:
            mask = [symbol != "H" for symbol in atoms.get_chemical_symbols()]
            atoms = atoms[mask]

        atoms.calc = self.calculator
        logger.info(
            "Атомов: %d, начальная энергия: %.2f eV",
            len(atoms),
            atoms.get_potential_energy(),
        )
        return atoms

    def minimize_energy(
        self,
        atoms: Atoms,
        fmax: float = 0.05,
        max_steps: int = 1000,
    ) -> Atoms:
        from ase.optimize import FIRE

        before = atoms.get_potential_energy()
        optimizer = FIRE(atoms, logfile=None)
        optimizer.run(fmax=fmax, steps=max_steps)
        after = atoms.get_potential_energy()
        logger.info("Минимизация: %.2f -> %.2f eV", before, after)
        return atoms

    def run_md(
        self,
        atoms: Atoms,
        temperature: float = 300.0,
        timestep: float = 1.0,
        num_steps: int = 10000,
        friction: float = 0.01,
        trajectory_file: str = "md_trajectory.dcd",
        log_interval: int = 100,
    ) -> dict:
        logger.info(
            "MD: T=%s K, dt=%s fs, шагов=%d",
            temperature,
            timestep,
            num_steps,
        )
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
        dyn = Langevin(
            atoms,
            timestep=timestep * units.fs,
            temperature_K=temperature,
            friction=friction,
        )

        energies = []
        temperatures = []
        trajectory = []

        def log_state():
            e_pot = atoms.get_potential_energy() / len(atoms)
            e_kin = atoms.get_kinetic_energy() / len(atoms)
            energies.append([e_pot, e_kin, e_pot + e_kin])
            temperatures.append(atoms.get_temperature())
            trajectory.append(atoms.get_positions().copy())

        dyn.attach(log_state, interval=log_interval)
        dyn.run(num_steps)

        self._save_trajectory(atoms, np.asarray(trajectory), trajectory_file)

        energies = np.asarray(energies)
        temperatures = np.asarray(temperatures)
        has_samples = energies.size > 0

        results = {
            "num_steps": num_steps,
            "energies": energies,
            "temperatures": temperatures,
            "trajectory_file": trajectory_file,
            "mean_temperature": float(np.mean(temperatures)) if has_samples else 0.0,
            "std_temperature": float(np.std(temperatures)) if has_samples else 0.0,
            "mean_energy": float(np.mean(energies[:, 2])) if has_samples else 0.0,
            "std_energy": float(np.std(energies[:, 2])) if has_samples else 0.0,
        }
        logger.info("Средняя температура: %.1f K", results["mean_temperature"])
        return results

    def _save_trajectory(self, atoms: Atoms, trajectory: np.ndarray, filename: str):
        path = Path(filename)

        if path.suffix == ".traj":
            from ase.io.trajectory import Trajectory

            traj = Trajectory(filename, "w")
            for positions in trajectory:
                frame = atoms.copy()
                frame.set_positions(positions)
                traj.write(frame)
            traj.close()
        elif path.suffix == ".xyz":
            with open(filename, "w", encoding="utf-8") as handle:
                for positions in trajectory:
                    frame = atoms.copy()
                    frame.set_positions(positions)
                    ase_write(handle, frame, format="xyz")
        elif path.suffix == ".pdb":
            with open(filename, "w", encoding="utf-8") as handle:
                for i, positions in enumerate(trajectory):
                    handle.write(f"MODEL     {i + 1:4d}\n")
                    frame = atoms.copy()
                    frame.set_positions(positions)
                    ase_write(handle, frame, format="proteindatabank")
                    handle.write("ENDMDL\n")
        elif path.suffix == ".dcd":
            self._save_dcd(atoms, trajectory, filename)
        else:
            raise ValueError(f"Неподдерживаемый формат траектории: {path.suffix}")

        logger.info("Траектория сохранена: %s (%d кадров)", filename, len(trajectory))

    @staticmethod
    def _save_dcd(atoms: Atoms, trajectory: np.ndarray, filename: str):
        import mdtraj as md
        import pandas as pd

        symbols = atoms.get_chemical_symbols()
        n_atoms = len(atoms)

        topology = md.Topology.from_dataframe(
            pd.DataFrame(
                {
                    "serial": np.arange(1, n_atoms + 1),
                    "name": symbols,
                    "element": symbols,
                    "resSeq": np.ones(n_atoms, dtype=int),
                    "resName": ["MOL"] * n_atoms,
                    "chainID": np.zeros(n_atoms, dtype=int),
                    "segmentID": [""] * n_atoms,
                }
            )
        )

        traj = md.Trajectory(
            xyz=np.asarray(trajectory, dtype=np.float64) / 10.0,
            topology=topology,
            time=np.arange(len(trajectory)),
        )
        traj.save_dcd(filename)
