"""
Простой докинг лиганда (геометрический + оптимизация).
Fallback если DiffDock недоступен.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
from typing import Tuple, List, Optional


class SimpleDocking:
    """
    Простой докинг:
    1. Генерирует 3D лиганда
    2. Находит центр масс белка (или активный сайт)
    3. Размещает лиганд с случайными ориентациями
    4. Оптимизирует геометрию
    """
    
    def __init__(self, num_poses: int = 10, seed: int = 42):
        self.num_poses = num_poses
        self.seed = seed
        np.random.seed(seed)
    
    def dock(
        self,
        protein_pdb: str,
        ligand_smiles: str,
        output_dir: str = "docking_results",
        binding_site_residues: Optional[List[int]] = None
    ) -> List[Tuple[str, float]]:
        """
        Выполняет докинг.
        
        Args:
            protein_pdb: Путь к PDB белка
            ligand_smiles: SMILES лиганда
            output_dir: Папка для результатов
            binding_site_residues: Номера остатков активного сайта (опц.)
        
        Returns:
            List[(путь_к_pdb, энергия)]
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        print(f"[Docking] Простой докинг...")
        print(f"          SMILES: {ligand_smiles}")
        
        # 1. Генерация 3D лиганда
        mol = self._prepare_ligand(ligand_smiles)
        
        # 2. Поиск центра связывания
        if binding_site_residues:
            center = self._get_binding_site_center(protein_pdb, binding_site_residues)
        else:
            center = self._get_protein_center(protein_pdb)
        
        print(f"          Центр: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
        
        # 3. Генерация поз
        poses = []
        
        for i in range(self.num_poses):
            # Клонируем молекулу
            mol_copy = Chem.Mol(mol)
            
            # Случайная ориентация
            mol_oriented = self._random_orientation(mol_copy)
            
            # Размещение в центре
            mol_placed = self._place_at_center(mol_oriented, center)
            
            # Оптимизация
            energy = self._optimize_geometry(mol_placed)
            
            # Сохранение
            pose_path = output_path / f"pose_{i+1}.pdb"
            Chem.MolToPDBFile(mol_placed, str(pose_path))
            
            poses.append((str(pose_path), energy))
        
        # Сортировка по энергии
        poses.sort(key=lambda x: x[1])
        
        print(f"[Docking] ✓ Сгенерировано {len(poses)} поз")
        print(f"          Лучшая энергия: {poses[0][1]:.2f} ккал/моль")
        
        return poses
    
    def _prepare_ligand(self, smiles: str) -> Chem.Mol:
        """Генерация 3D конформера."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Невалидный SMILES: {smiles}")
        
        # Добавление водородов
        mol = Chem.AddHs(mol)
        
        # Генерация 3D
        params = AllChem.ETKDGv3()
        params.randomSeed = self.seed

        if AllChem.EmbedMolecule(mol, params) != 0:
            # Fallback
            if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
                raise ValueError(
                    f"Не удалось сгенерировать 3D-структуру для SMILES: {smiles}"
                )

        # Оптимизация (если MMFF недоступен - пробуем UFF, иначе пропускаем)
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                pass

        return mol
    
    @staticmethod
    def _get_protein_center(pdb_path: str) -> np.ndarray:
        """Центр масс всех атомов белка."""
        coords = []
        
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
        
        if not coords:
            raise ValueError("Нет атомов в PDB файле")
        
        return np.mean(coords, axis=0)
    
    @staticmethod
    def _get_binding_site_center(pdb_path: str, residues: List[int]) -> np.ndarray:
        """Центр активного сайта."""
        coords = []
        
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM"):
                    resnum = int(line[22:26].strip())
                    if resnum in residues:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
        
        if not coords:
            raise ValueError(f"Остатки {residues} не найдены в PDB")
        
        return np.mean(coords, axis=0)
    
    def _random_orientation(self, mol: Chem.Mol) -> Chem.Mol:
        """Случайная ориентация молекулы."""
        conf = mol.GetConformer()
        
        # Случайная матрица поворота
        # Используем случайные углы Эйлера
        alpha = np.random.uniform(0, 2*np.pi)
        beta = np.random.uniform(0, 2*np.pi)
        gamma = np.random.uniform(0, 2*np.pi)
        
        # Матрица поворота
        Rz_alpha = np.array([
            [np.cos(alpha), -np.sin(alpha), 0],
            [np.sin(alpha), np.cos(alpha), 0],
            [0, 0, 1]
        ])
        
        Ry_beta = np.array([
            [np.cos(beta), 0, np.sin(beta)],
            [0, 1, 0],
            [-np.sin(beta), 0, np.cos(beta)]
        ])
        
        Rz_gamma = np.array([
            [np.cos(gamma), -np.sin(gamma), 0],
            [np.sin(gamma), np.cos(gamma), 0],
            [0, 0, 1]
        ])
        
        R = Rz_gamma @ Ry_beta @ Rz_alpha
        
        # Применение поворота
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            pos_vec = np.array([pos.x, pos.y, pos.z])
            new_pos = R @ pos_vec
            conf.SetAtomPosition(i, tuple(new_pos))
        
        return mol
    
    def _place_at_center(self, mol: Chem.Mol, center: np.ndarray) -> Chem.Mol:
        """Размещает молекулу в заданном центре."""
        conf = mol.GetConformer()
        
        # Текущий центр лиганда
        lig_coords = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            lig_coords.append([pos.x, pos.y, pos.z])
        
        lig_center = np.mean(lig_coords, axis=0)
        
        # Сдвиг
        shift = center - lig_center
        
        # Добавляем небольшую случайность (±2 Å)
        shift += np.random.randn(3) * 2.0
        
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            new_pos = [
                pos.x + shift[0],
                pos.y + shift[1],
                pos.z + shift[2]
            ]
            conf.SetAtomPosition(i, tuple(new_pos))
        
        return mol
    
    @staticmethod
    def _optimize_geometry(mol: Chem.Mol, max_iters: int = 500) -> float:
        """Оптимизация геометрии (MMFF)."""
        props = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, props) if props is not None else None

        if ff is None:
            # Fallback на UFF
            ff = AllChem.UFFGetMoleculeForceField(mol)

        if ff is None:
            # Нет доступного силового поля - возвращаем нулевую энергию
            return 0.0

        ff.Minimize(maxIts=max_iters)
        energy = ff.CalcEnergy()

        return energy
    
    def create_complex(
        self,
        protein_pdb: str,
        ligand_pdb: str,
        output_pdb: str
    ):
        """Объединяет белок и лиганд в один PDB."""
        with open(output_pdb, "w") as out:
            # Белок
            with open(protein_pdb) as f:
                for line in f:
                    if line.startswith("ATOM"):
                        out.write(line)
            
            out.write("TER\n")
            
            # Лиганд (меняем цепь на L)
            with open(ligand_pdb) as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        # Замена цепи
                        line = line[:21] + "L" + line[22:]
                        # Замена остатка
                        line = line[:17] + "LIG" + line[20:]
                        out.write(line)
            
            out.write("END\n")
        
        print(f"[Docking] ✓ Комплекс: {output_pdb}")


# ============ ТЕСТ ============
if __name__ == "__main__":
    docker = SimpleDocking(num_poses=5)
    
    poses = docker.dock(
        protein_pdb="test_data/protein.pdb",
        ligand_smiles="CC(=O)OC1=CC=CC=C1C(=O)O",
        output_dir="test_docking"
    )
    
    # Создание комплекса с лучшей позой
    docker.create_complex(
        "test_data/protein.pdb",
        poses[0][0],
        "test_docking/complex.pdb"
    )
