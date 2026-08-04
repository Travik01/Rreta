"""
MACE-based молекулярная динамика.
Использует предобученный ML-потенциал для высокоточной симуляции.
"""

import numpy as np
from pathlib import Path
from typing import Optional, Dict, Tuple
from ase import Atoms
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase import units
from ase.io import read as ase_read
from ase.io import write as ase_write

# torch и mace-torch загружаются лениво (в __init__),
# чтобы модуль импортировался и без установленных тяжёлых зависимостей.
try:
    import torch
    _TORCH_AVAILABLE = True
except ImportError:
    torch = None
    _TORCH_AVAILABLE = False

try:
    from mace.calculators import MACECalculator
    _MACE_AVAILABLE = True
except ImportError:
    MACECalculator = None
    _MACE_AVAILABLE = False


class MACESimulator:
    """
    MD-симулятор на базе MACE ML-потенциала.
    
    Преимущества:
    - DFT-уровень точности (~1 ккал/моль)
    - Быстрее DFT в 1000+ раз
    - Предобучена на Materials Project (1.5M структур)
    - Переносится на биомолекулы
    """
    
    def __init__(
        self,
        model_path: str = "mace_mp_0",  # Предобученная модель
        device: str = "cuda",
        default_dtype: str = "float32"
    ):
        if not _MACE_AVAILABLE or not _TORCH_AVAILABLE:
            raise ImportError(
                "Для MACE-симуляции установите зависимости:\n"
                "    pip install mace-torch\n"
                "(потребуется PyTorch с поддержкой CUDA для GPU)"
            )

        if device == "cuda" and not torch.cuda.is_available():
            device = "cpu"
        self.device = device
        
        print(f"[MACE] Инициализация...")
        print(f"       Device: {self.device}")
        
        # Загрузка MACE калькулятора
        try:
            self.calculator = MACECalculator(
                model_paths=model_path,
                device=self.device,
                default_dtype=default_dtype
            )
            print(f"[MACE] ✓ Модель загружена: {model_path}")
        
        except Exception as e:
            print(f"[MACE] ⚠ Ошибка загрузки модели: {e}")
            print(f"[MACE]   Попытка загрузки из HuggingFace...")
            
            # Fallback: загрузка из интернета
            self.calculator = MACECalculator(
                model_paths="https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0/mace_mp_0.model",
                device=self.device,
                default_dtype=default_dtype
            )
            print(f"[MACE] ✓ Модель загружена из репозитория")
    
    def prepare_system(
        self,
        pdb_file: str,
        remove_hydrogens: bool = False
    ) -> Atoms:
        """
        Подготавливает систему из PDB.
        
        Args:
            pdb_file: Путь к PDB файлу
            remove_hydrogens: Удалить водороды (MACE работает без них)
        
        Returns:
            ASE Atoms объект
        """
        print(f"[MACE] Загрузка структуры: {pdb_file}")
        
        # Загрузка через ASE
        atoms = ase_read(pdb_file)
        
        # Удаление водородов (опционально)
        if remove_hydrogens:
            mask = atoms.get_chemical_symbols() != 'H'
            atoms = atoms[mask]
            print(f"       Водороды удалены")
        
        # Установка калькулятора
        atoms.calc = self.calculator
        
        print(f"       Атомов: {len(atoms)}")
        print(f"       Формула: {atoms.get_chemical_formula()}")
        
        # Проверка энергии
        energy = atoms.get_potential_energy()
        print(f"       Начальная энергия: {energy:.2f} eV")
        
        return atoms
    
    def minimize_energy(
        self,
        atoms: Atoms,
        fmax: float = 0.05,  # eV/Å
        max_steps: int = 1000
    ) -> Atoms:
        """
        Минимизация энергии.
        
        Args:
            atoms: ASE Atoms
            fmax: Максимальная сила (критерий сходимости)
            max_steps: Максимум шагов
        
        Returns:
            Оптимизированная структура
        """
        from ase.optimize import FIRE
        
        print(f"[MACE] Минимизация энергии (fmax={fmax} eV/Å)...")
        
        energy_before = atoms.get_potential_energy()
        
        optimizer = FIRE(atoms, logfile=None)
        optimizer.run(fmax=fmax, steps=max_steps)
        
        energy_after = atoms.get_potential_energy()
        
        print(f"       Энергия до:  {energy_before:.2f} eV")
        print(f"       Энергия после: {energy_after:.2f} eV")
        print(f"       Δ Энергия: {energy_after - energy_before:.2f} eV")
        
        return atoms
    
    def run_md(
        self,
        atoms: Atoms,
        temperature: float = 300.0,  # K
        timestep: float = 1.0,       # fs
        num_steps: int = 10000,
        friction: float = 0.01,      # 1/fs (для Langevin)
        trajectory_file: str = "md_trajectory.dcd",
        log_interval: int = 100
    ) -> Dict:
        """
        Запуск MD-симуляции.
        
        Args:
            atoms: ASE Atoms
            temperature: Температура (K)
            timestep: Шаг времени (fs)
            num_steps: Количество шагов
            friction: Коэффициент трения для термостата Langevin
            trajectory_file: Файл траектории
            log_interval: Интервал логирования
        
        Returns:
            Словарь с результатами
        """
        print(f"[MACE] Запуск MD-симуляции...")
        print(f"       T = {temperature} K")
        print(f"       dt = {timestep} fs")
        print(f"       Шагов = {num_steps} ({num_steps * timestep / 1000:.2f} пс)")
        
        # Инициализация скоростей (Maxwell-Boltzmann)
        MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
        
        # Langevin термостат (NVT)
        dyn = Langevin(
            atoms,
            timestep=timestep * units.fs,
            temperature_K=temperature,
            friction=friction
        )
        
        # Логирование
        energies = []
        temperatures = []
        trajectory = []
        
        def log_state():
            e_pot = atoms.get_potential_energy() / len(atoms)  # eV/atom
            e_kin = atoms.get_kinetic_energy() / len(atoms)
            e_tot = e_pot + e_kin
            temp = atoms.get_temperature()
            
            energies.append([e_pot, e_kin, e_tot])
            temperatures.append(temp)
            trajectory.append(atoms.get_positions().copy())
        
        # Attach логирования
        dyn.attach(log_state, interval=log_interval)
        
        # Запуск
        print(f"[MACE] Симуляция...")
        dyn.run(num_steps)
        
        print(f"[MACE] ✓ Симуляция завершена")
        
        # Сохранение траектории
        self._save_trajectory(
            atoms,
            np.array(trajectory),
            trajectory_file
        )
        
        # Статистика
        energies = np.array(energies)
        temperatures = np.array(temperatures)
        
        results = {
            "num_steps": num_steps,
            "energies": energies,
            "temperatures": temperatures,
            "trajectory_file": trajectory_file,
            "mean_temperature": np.mean(temperatures),
            "std_temperature": np.std(temperatures),
            "mean_energy": np.mean(energies[:, 2]),
            "std_energy": np.std(energies[:, 2])
        }
        
        print(f"       <T> = {results['mean_temperature']:.1f} ± {results['std_temperature']:.1f} K")
        print(f"       <E> = {results['mean_energy']:.4f} ± {results['std_energy']:.4f} eV/atom")
        
        return results
    
    def _save_trajectory(
        self,
        atoms: Atoms,
        trajectory: np.ndarray,
        filename: str
    ):
        """Сохранение траектории в различных форматах."""
        path = Path(filename)
        
        # ASE trajectory format
        if path.suffix == ".traj":
            from ase.io.trajectory import Trajectory
            
            traj = Trajectory(filename, 'w')
            for positions in trajectory:
                atoms_copy = atoms.copy()
                atoms_copy.set_positions(positions)
                traj.write(atoms_copy)
            traj.close()
        
        # XYZ format
        elif path.suffix == ".xyz":
            with open(filename, "w") as f:
                for positions in trajectory:
                    atoms_copy = atoms.copy()
                    atoms_copy.set_positions(positions)
                    ase_write(f, atoms_copy, format="xyz")
        
        # PDB multi-model
        elif path.suffix == ".pdb":
            with open(filename, "w") as f:
                for i, positions in enumerate(trajectory):
                    f.write(f"MODEL     {i+1:4d}\n")
                    atoms_copy = atoms.copy()
                    atoms_copy.set_positions(positions)
                    ase_write(f, atoms_copy, format="proteindatabank")
                    f.write("ENDMDL\n")

        # DCD (рекомендуется: читается MDTraj/анализатором)
        elif path.suffix == ".dcd":
            self._save_dcd(atoms, trajectory, filename)

        else:
            raise ValueError(
                f"Неизвестный формат траектории: {path.suffix} "
                f"(поддерживаются: .traj, .xyz, .pdb, .dcd)"
            )

        print(f"[MACE] ✓ Траектория сохранена: {filename} ({len(trajectory)} кадров)")

    @staticmethod
    def _save_dcd(atoms: Atoms, trajectory: np.ndarray, filename: str):
        """
        Сохранение траектории в формате DCD через MDTraj
        (единицы: ASE/Å -> MDTraj/нм).
        """
        import pandas as pd
        import mdtraj as md

        symbols = atoms.get_chemical_symbols()
        n_atoms = len(atoms)

        # Топология из ASE-структуры (одинаковый порядок атомов,
        # как в исходном PDB-файле)
        df = pd.DataFrame({
            "serial": np.arange(1, n_atoms + 1),
            "name": symbols,
            "element": symbols,
            "resSeq": np.ones(n_atoms, dtype=int),
            "resName": ["MOL"] * n_atoms,
            "chainID": np.zeros(n_atoms, dtype=int),
            "segmentID": [""] * n_atoms,
        })
        topology = md.Topology.from_dataframe(df)

        xyz_nm = np.asarray(trajectory, dtype=np.float64) / 10.0  # Å -> нм
        traj = md.Trajectory(
            xyz=xyz_nm,
            topology=topology,
            time=np.arange(len(trajectory))
        )
        traj.save_dcd(filename)


# ============ ТЕСТ ============
if __name__ == "__main__":
    # Инициализация
    simulator = MACESimulator(device="cuda")
    
    # Загрузка структуры (тестовая молекула)
    # Можно использовать свой PDB файл
    from ase.build import molecule
    atoms = molecule('CH3CH2OH')  # Этанол для теста
    atoms.calc = simulator.calculator
    
    # Минимизация
    atoms = simulator.minimize_energy(atoms, fmax=0.01)
    
    # MD-симуляция
    results = simulator.run_md(
        atoms,
        temperature=300.0,
        timestep=0.5,
        num_steps=5000,
        trajectory_file="test_mace_md.traj"
    )
    
    print("\nРезультаты:")
    for k, v in results.items():
        if not isinstance(v, np.ndarray):
            print(f"  {k}: {v}")
