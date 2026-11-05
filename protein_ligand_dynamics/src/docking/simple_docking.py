"""Геометрический докинг лиганда в структуру белка."""

import logging
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)


class SimpleDocking:
    def __init__(self, num_poses: int = 10, seed: int = 42):
        self.num_poses = num_poses
        self.rng = np.random.default_rng(seed)

    def dock(
        self,
        protein_pdb: str,
        ligand_smiles: str,
        output_dir: str = "docking_results",
        binding_site_residues: list[int] | None = None,
    ) -> list[tuple[str, float]]:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        logger.info("Докинг %s, поз: %d", ligand_smiles, self.num_poses)
        mol = self._prepare_ligand(ligand_smiles)

        if binding_site_residues:
            center = self._get_binding_site_center(protein_pdb, binding_site_residues)
        else:
            center = self._get_protein_center(protein_pdb)

        poses = []
        for i in range(self.num_poses):
            oriented = self._random_orientation(Chem.Mol(mol))
            placed = self._place_at_center(oriented, center)
            energy = self._optimize_geometry(placed)
            pose_path = output_path / f"pose_{i + 1}.pdb"
            Chem.MolToPDBFile(placed, str(pose_path))
            poses.append((str(pose_path), energy))

        poses.sort(key=lambda item: item[1])
        if poses:
            logger.info("Поз: %d, лучшая энергия: %.2f", len(poses), poses[0][1])

        return poses

    def _prepare_ligand(self, smiles: str) -> Chem.Mol:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Некорректный SMILES: {smiles}")

        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = int(self.rng.integers(0, 2**31 - 1))

        if (
            AllChem.EmbedMolecule(mol, params) != 0
            and AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0
        ):
            raise ValueError(f"Не удалось построить 3D-структуру для SMILES: {smiles}")

        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except (RuntimeError, ValueError):
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            except (RuntimeError, ValueError):
                logger.warning("Оптимизация конформера недоступна для %s", smiles)

        return mol

    @staticmethod
    def _get_protein_center(pdb_path: str) -> np.ndarray:
        coords = []
        with open(pdb_path, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("ATOM"):
                    coords.append(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    )

        if not coords:
            raise ValueError("В PDB нет атомов белка")

        return np.mean(coords, axis=0)

    @staticmethod
    def _get_binding_site_center(pdb_path: str, residues: list[int]) -> np.ndarray:
        coords = []
        with open(pdb_path, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("ATOM") and int(line[22:26].strip()) in residues:
                    coords.append(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    )

        if not coords:
            raise ValueError(f"Остатки {residues} не найдены в PDB")

        return np.mean(coords, axis=0)

    def _random_orientation(self, mol: Chem.Mol) -> Chem.Mol:
        alpha, beta, gamma = self.rng.uniform(0, 2 * np.pi, 3)

        rz_alpha = np.array(
            [[np.cos(alpha), -np.sin(alpha), 0], [np.sin(alpha), np.cos(alpha), 0], [0, 0, 1]]
        )
        ry_beta = np.array(
            [[np.cos(beta), 0, np.sin(beta)], [0, 1, 0], [-np.sin(beta), 0, np.cos(beta)]]
        )
        rz_gamma = np.array(
            [[np.cos(gamma), -np.sin(gamma), 0], [np.sin(gamma), np.cos(gamma), 0], [0, 0, 1]]
        )
        rotation = rz_gamma @ ry_beta @ rz_alpha

        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            new_pos = rotation @ np.array([pos.x, pos.y, pos.z])
            conf.SetAtomPosition(i, tuple(new_pos))

        return mol

    def _place_at_center(self, mol: Chem.Mol, center: np.ndarray) -> Chem.Mol:
        conf = mol.GetConformer()
        coords = np.array(
            [[p.x, p.y, p.z] for p in (conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms()))]
        )
        shift = center - coords.mean(axis=0) + self.rng.normal(0, 2.0, 3)

        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pos.x + shift[0], pos.y + shift[1], pos.z + shift[2]))

        return mol

    @staticmethod
    def _optimize_geometry(mol: Chem.Mol, max_iters: int = 500) -> float:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        force_field = AllChem.MMFFGetMoleculeForceField(mol, props) if props is not None else None

        if force_field is None:
            force_field = AllChem.UFFGetMoleculeForceField(mol)
        if force_field is None:
            return 0.0

        force_field.Minimize(maxIts=max_iters)
        return float(force_field.CalcEnergy())

    def create_complex(self, protein_pdb: str, ligand_pdb: str, output_pdb: str):
        with open(output_pdb, "w", encoding="utf-8") as out:
            with open(protein_pdb, encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("ATOM"):
                        out.write(line)

            out.write("TER\n")

            with open(ligand_pdb, encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith(("ATOM", "HETATM")):
                        line = line[:21] + "L" + line[22:]
                        line = line[:17] + "LIG" + line[20:]
                        out.write(line)

            out.write("END\n")

        logger.info("Комплекс сохранён: %s", output_pdb)
