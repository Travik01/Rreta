"""
Тесты для модуля докинга.
"""

import pytest
import tempfile
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from docking.simple_docking import SimpleDocking


class TestSimpleDocking:
    """Тесты для простого докинга."""
    
    @pytest.fixture
    def docker(self):
        """Фикстура: докер с фиксированным seed."""
        return SimpleDocking(num_poses=3, seed=42)
    
    @pytest.fixture
    def test_protein(self, tmp_path):
        """Создание тестового PDB белка."""
        pdb_content = """
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM      3  C   ALA A   1      12.000  12.000  12.000  1.00  0.00           C
ATOM      4  N   GLY A   2      13.000  13.000  13.000  1.00  0.00           N
ATOM      5  CA  GLY A   2      14.000  14.000  14.000  1.00  0.00           C
"""
        pdb_file = tmp_path / "protein.pdb"
        pdb_file.write_text(pdb_content)
        return str(pdb_file)
    
    def test_ligand_preparation(self, docker):
        """Тест генерации 3D лиганда из SMILES."""
        smiles = "CC(=O)O"  # Уксусная кислота
        mol = docker._prepare_ligand(smiles)
        
        assert mol is not None
        assert mol.GetNumAtoms() > 0
    
    def test_protein_center_calculation(self, docker, test_protein):
        """Тест расчёта центра белка."""
        center = docker._get_protein_center(test_protein)
        
        assert len(center) == 3
        assert center[0] == pytest.approx(12.0, abs=0.1)
        assert center[1] == pytest.approx(12.0, abs=0.1)
        assert center[2] == pytest.approx(12.0, abs=0.1)
    
    def test_docking_workflow(self, docker, test_protein):
        """Тест полного процесса докинга."""
        with tempfile.TemporaryDirectory() as tmpdir:
            poses = docker.dock(
                protein_pdb=test_protein,
                ligand_smiles="CC(=O)O",
                output_dir=tmpdir
            )
            
            # Проверки
            assert len(poses) == 3
            
            for pose_file, energy in poses:
                assert Path(pose_file).exists()
                assert isinstance(energy, float)
            
            # Энергии должны быть отсортированы
            energies = [e for _, e in poses]
            assert energies == sorted(energies)
    
    def test_complex_creation(self, docker, test_protein, tmp_path):
        """Тест создания комплекса."""
        # Создаём тестовый лиганд PDB
        ligand_pdb = tmp_path / "ligand.pdb"
        ligand_content = """
HETATM    1  C   LIG L   1      15.000  15.000  15.000  1.00  0.00           C
HETATM    2  O   LIG L   1      16.000  16.000  16.000  1.00  0.00           O
"""
        ligand_pdb.write_text(ligand_content)
        
        complex_pdb = tmp_path / "complex.pdb"
        
        docker.create_complex(
            test_protein,
            str(ligand_pdb),
            str(complex_pdb)
        )
        
        # Проверки
        assert complex_pdb.exists()
        
        content = complex_pdb.read_text()
        assert "ATOM" in content  # Белок
        assert "TER" in content
        assert "HETATM" in content or "ATOM" in content  # Лиганд
        assert "LIG" in content  # Имя лиганда


# ============ ЗАПУСК ============

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
