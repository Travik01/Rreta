"""
Анализ MD-траекторий: RMSD, RMSF, Rg, контакты, энергия связывания.
"""

import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from scipy.spatial.distance import cdist
import pandas as pd


class TrajectoryAnalyzer:
    """
    Комплексный анализ MD-траектории.
    
    Метрики:
    - RMSD (стабильность структуры)
    - RMSF (флуктуации атомов)
    - Radius of Gyration (компактность)
    - Межмолекулярные контакты
    - Водородные связи
    - Расстояния белок-лиганд
    """
    
    def __init__(self, trajectory_file: str, topology_file: str):
        """
        Args:
            trajectory_file: Путь к траектории (.dcd, .xtc, .traj)
            topology_file: Путь к топологии (.pdb)
        """
        print(f"[Analyzer] Загрузка траектории...")
        
        self.traj = md.load(trajectory_file, top=topology_file)
        self.topology = self.traj.topology
        
        print(f"           Кадров: {self.traj.n_frames}")
        print(f"           Атомов: {self.traj.n_atoms}")
        print(f"           Время: {self.traj.time[-1]:.2f} пс")
        
        # Выделение белка и лиганда
        self.protein_indices = self.topology.select("protein")
        self.ligand_indices = self.topology.select("resname LIG or resname UNK")
        
        if len(self.ligand_indices) == 0:
            # Попытка найти лиганд по цепи
            self.ligand_indices = self.topology.select("chainid 1")
        
        print(f"           Атомов белка: {len(self.protein_indices)}")
        print(f"           Атомов лиганда: {len(self.ligand_indices)}")
    
    def compute_rmsd(
        self,
        reference_frame: int = 0,
        selection: str = "protein and name CA",
        align: bool = True
    ) -> np.ndarray:
        """
        Вычисляет RMSD относительно референсной структуры.
        
        Args:
            reference_frame: Номер референсного кадра
            selection: Атомы для расчёта (MDTraj DSL)
            align: Выравнивать структуры перед расчётом
        
        Returns:
            RMSD для каждого кадра (в Ангстремах)
        """
        print(f"[RMSD] Расчёт...")
        
        atom_indices = self.topology.select(selection)
        
        if len(atom_indices) == 0:
            raise ValueError(f"Нет атомов для селекции: {selection}")
        
        # Выравнивание
        if align:
            self.traj.superpose(self.traj, frame=reference_frame, atom_indices=atom_indices)
        
        # RMSD
        rmsd = md.rmsd(self.traj, self.traj, frame=reference_frame, atom_indices=atom_indices)
        rmsd *= 10  # нм → Å
        
        print(f"       RMSD mean: {np.mean(rmsd):.2f} ± {np.std(rmsd):.2f} Å")
        print(f"       RMSD max: {np.max(rmsd):.2f} Å")
        
        return rmsd
    
    def compute_rmsf(
        self,
        selection: str = "protein and name CA"
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Root Mean Square Fluctuation - флуктуации каждого атома.
        
        Returns:
            (rmsf_values, residue_indices)
        """
        print(f"[RMSF] Расчёт...")
        
        atom_indices = self.topology.select(selection)
        
        # RMSF: флуктуации атомов относительно среднего положения
        # (x[i] = sqrt( mean_t || r_i(t) - <r_i> ||^2 ) ), нм -> Å
        xyz = self.traj.xyz[:, atom_indices]                      # (frames, atoms, 3)
        mean_pos = xyz.mean(axis=0, keepdims=True)                # средняя структура
        rmsf = np.sqrt(np.mean(np.sum((xyz - mean_pos) ** 2, axis=2), axis=0)) * 10
        
        # Получение номеров остатков
        residue_indices = [self.topology.atom(i).residue.resSeq for i in atom_indices]
        
        print(f"       RMSF mean: {np.mean(rmsf):.2f} Å")
        print(f"       RMSF max: {np.max(rmsf):.2f} Å (остаток {residue_indices[np.argmax(rmsf)]})")
        
        return rmsf, np.array(residue_indices)
    
    def compute_radius_of_gyration(
        self,
        selection: str = "protein"
    ) -> np.ndarray:
        """
        Радиус гирации - мера компактности белка.
        
        Returns:
            Rg для каждого кадра (в Ангстремах)
        """
        print(f"[Rg] Расчёт радиуса гирации...")
        
        atom_indices = self.topology.select(selection)
        # Важно: atom_slice() выбирает атомы, а traj[индексы] - кадры!
        rg = md.compute_rg(self.traj.atom_slice(atom_indices)) * 10  # нм → Å
        
        print(f"     Rg mean: {np.mean(rg):.2f} ± {np.std(rg):.2f} Å")
        
        return rg
    
    def compute_contacts(
        self,
        cutoff: float = 0.45,  # нм (4.5 Å)
        scheme: str = "closest-heavy"  # сохранено для совместимости API
    ) -> Tuple[np.ndarray, List]:
        """
        Межмолекулярные контакты белок-лиганд.

        Контакт = атом белка на расстоянии < cutoff от любого атома лиганда.
        Считается напрямую по атомным расстояниям (надёжнее, чем
        md.compute_contacts, который пропускает нестандартные остатки).

        Returns:
            (contact_frequencies, contact_pairs) - частоты по остаткам белка
        """
        print(f"[Contacts] Расчёт контактов (cutoff={cutoff*10:.1f} Å)...")

        if len(self.ligand_indices) == 0:
            print("       ⚠ Лиганд не найден, пропуск")
            return np.array([]), []

        # Атомы белка, сгруппированные по остаткам
        protein_atoms = list(self.protein_indices)
        residue_of_atom = {i: self.topology.atom(i).residue for i in protein_atoms}
        residues_in_order = []
        for res in self.topology.residues:
            if res.index in {r.index for r in residue_of_atom.values()}:
                residues_in_order.append(res)

        # Минимальное расстояние белок-лиганд по кадрам
        n_frames = len(self.traj)
        min_dist = np.full((n_frames, len(protein_atoms)), np.inf)

        for f, frame in enumerate(self.traj):
            prot_coords = frame.xyz[0, self.protein_indices]
            lig_coords = frame.xyz[0, self.ligand_indices]
            dist_matrix = cdist(prot_coords, lig_coords)  # нм
            min_dist[f] = dist_matrix.min(axis=1)

        # Частота контакта для каждого остатка белка
        contact_pairs = []
        contact_frequencies = []
        for res in residues_in_order:
            atom_positions = [i for i in protein_atoms if residue_of_atom[i] is res]
            res_min = min_dist[:, atom_positions]  # (frames, atoms_res)
            freq = float(np.mean(np.any(res_min < cutoff, axis=1)))
            if freq > 0:
                contact_pairs.append((res.index, "LIG"))
                contact_frequencies.append(freq)

        contact_frequencies = np.array(contact_frequencies) if contact_frequencies else np.array([])

        if len(contact_frequencies) == 0:
            print("       ⚠ Нет межмолекулярных контактов")
            return np.array([]), []

        # Топ-5 контактов
        top_indices = np.argsort(contact_frequencies)[-5:][::-1]

        print(f"       Найдено {len(contact_frequencies)} остатков в контакте с лигандом")
        print(f"       Топ-5 контактов:")

        for idx in top_indices:
            r1 = contact_pairs[idx][0]
            res1 = self.topology.residue(r1)
            freq = contact_frequencies[idx]
            print(f"         {res1.name}{res1.resSeq} - LIG: {freq*100:.1f}%")

        return contact_frequencies, contact_pairs
    
    def compute_ligand_protein_distance(self) -> np.ndarray:
        """
        Минимальное расстояние между лигандом и белком.
        
        Returns:
            Расстояние для каждого кадра (Å)
        """
        print(f"[Distance] Расчёт расстояния белок-лиганд...")
        
        if len(self.ligand_indices) == 0:
            raise ValueError("Лиганд не найден")
        
        distances = []
        
        for frame in self.traj:
            protein_coords = frame.xyz[0, self.protein_indices]
            ligand_coords = frame.xyz[0, self.ligand_indices]
            
            # Матрица расстояний
            dist_matrix = cdist(protein_coords, ligand_coords)
            min_dist = np.min(dist_matrix)
            
            distances.append(min_dist * 10)  # нм → Å
        
        distances = np.array(distances)
        
        print(f"       Расстояние mean: {np.mean(distances):.2f} ± {np.std(distances):.2f} Å")
        print(f"       Расстояние min: {np.min(distances):.2f} Å")
        
        return distances
    
    def compute_hydrogen_bonds(
        self,
        distance_cutoff: float = 0.25,  # нм
        angle_cutoff: float = 120.0      # градусы
    ) -> Tuple[List, np.ndarray]:
        """
        Водородные связи между белком и лигандом.
        
        Returns:
            (hbond_pairs, hbond_frequencies)
        """
        print(f"[H-bonds] Расчёт водородных связей...")
        
        if len(self.ligand_indices) == 0:
            return [], np.array([])
        
        # Вычисление H-связей
        hbonds = md.baker_hubbard(
            self.traj,
            distance_cutoff=distance_cutoff,
            angle_cutoff=angle_cutoff,
            periodic=False
        )
        
        if len(hbonds) == 0:
            print("       Водородные связи не найдены")
            return [], np.array([])
        
        # Фильтрация: только межмолекулярные
        protein_set = set(self.protein_indices)
        ligand_set = set(self.ligand_indices)
        
        inter_hbonds = []
        
        for donor, h, acceptor in hbonds:
            if (donor in protein_set and acceptor in ligand_set) or \
               (donor in ligand_set and acceptor in protein_set):
                inter_hbonds.append((donor, h, acceptor))
        
        if not inter_hbonds:
            print("       Межмолекулярные H-связи не найдены")
            return [], np.array([])
        
        # Частота H-связей (упрощённо - считаем уникальные пары)
        unique_pairs = list(set([(d, a) for d, h, a in inter_hbonds]))
        frequencies = np.ones(len(unique_pairs))  # Заглушка
        
        print(f"       Найдено {len(unique_pairs)} уникальных H-связей")
        
        return unique_pairs, frequencies
    
    def generate_report(self, output_dir: str = "analysis_results") -> Dict:
        """
        Генерирует полный отчёт с графиками.
        
        Returns:
            Словарь со всеми метриками
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        print(f"\n{'='*60}")
        print(f"  ПОЛНЫЙ АНАЛИЗ ТРАЕКТОРИИ")
        print(f"{'='*60}\n")
        
        results = {}
        
        # 1. RMSD
        rmsd = self.compute_rmsd()
        results["rmsd"] = rmsd
        self._plot_timeseries(
            rmsd,
            "RMSD (Å)",
            "RMSD: Стабильность структуры",
            output_path / "rmsd.png"
        )
        
        # 2. RMSF
        rmsf, residues = self.compute_rmsf()
        results["rmsf"] = rmsf
        results["rmsf_residues"] = residues
        self._plot_per_residue(
            residues,
            rmsf,
            "RMSF (Å)",
            "RMSF: Флуктуации остатков",
            output_path / "rmsf.png"
        )
        
        # 3. Radius of Gyration
        rg = self.compute_radius_of_gyration()
        results["radius_of_gyration"] = rg
        self._plot_timeseries(
            rg,
            "Rg (Å)",
            "Радиус гирации",
            output_path / "rg.png"
        )
        
        # 4. Контакты
        if len(self.ligand_indices) > 0:
            contact_freq, contact_pairs = self.compute_contacts()
            results["contact_frequencies"] = contact_freq
            results["contact_pairs"] = contact_pairs
            
            # 5. Расстояние белок-лиганд
            distances = self.compute_ligand_protein_distance()
            results["ligand_distance"] = distances
            self._plot_timeseries(
                distances,
                "Расстояние (Å)",
                "Белок-лиганд расстояние",
                output_path / "distance.png"
            )
            
            # 6. H-связи
            hbonds, hbond_freq = self.compute_hydrogen_bonds()
            results["hydrogen_bonds"] = hbonds
        
        # Сохранение CSV
        self._save_csv(results, output_path / "metrics.csv")
        
        print(f"\n{'='*60}")
        print(f"  ✓ Анализ завершён")
        print(f"  Результаты: {output_path}/")
        print(f"{'='*60}\n")
        
        return results
    
    @staticmethod
    def _plot_timeseries(data, ylabel, title, filename):
        """Построение временного ряда."""
        plt.figure(figsize=(10, 4))
        plt.plot(data, linewidth=0.8)
        plt.xlabel("Кадр")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
    
    @staticmethod
    def _plot_per_residue(residues, data, ylabel, title, filename):
        """График по остаткам."""
        plt.figure(figsize=(12, 4))
        plt.bar(residues, data, width=1.0, alpha=0.7)
        plt.xlabel("Номер остатка")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
    
    @staticmethod
    def _save_csv(results, filename):
        """Сохранение метрик в CSV (массивы разной длины выравниваются NaN)."""
        data = {}
        
        for key, value in results.items():
            if isinstance(value, np.ndarray) and value.ndim == 1 and len(value) > 0:
                data[key] = np.asarray(value, dtype=float)
        
        if data:
            # Выравнивание по максимальной длине
            max_len = max(len(v) for v in data.values())
            aligned = {
                key: np.pad(v, (0, max_len - len(v)), constant_values=np.nan)
                for key, v in data.items()
            }
            df = pd.DataFrame(aligned)
            df.to_csv(filename, index=False)
            print(f"       ✓ Метрики сохранены: {filename}")


# ============ ТЕСТ ============
if __name__ == "__main__":
    analyzer = TrajectoryAnalyzer(
        trajectory_file="test_data/trajectory.dcd",
        topology_file="test_data/complex.pdb"
    )
    
    results = analyzer.generate_report(output_dir="test_analysis")
