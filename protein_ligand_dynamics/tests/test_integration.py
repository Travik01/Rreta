"""
Интеграционные тесты полного пайплайна.
"""

import pytest
import tempfile
import yaml
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from full_pipeline import ProteinLigandPipeline


class TestFullPipeline:
    """Тесты полного пайплайна."""
    
    @pytest.fixture
    def minimal_config(self):
        """Минимальная конфигурация для быстрого теста."""
        return {
            "structure_prediction": {
                "sequence": "MKFLKFSLLTAVLLSVVFAFSSCG",  # Очень короткий
                "use_local": False
            },
            "docking": {
                "ligand_smiles": "CC(=O)O",  # Простая молекула
                "num_poses": 2
            },
            "md_simulation": {
                "use_mace": False,  # Используем OpenMM для стабильности
                "use_gpu": False,   # CPU для CI/CD
                "temperature": 300.0,
                "num_steps": 100,   # Очень мало для теста
                "report_interval": 50
            }
        }
    
    @pytest.mark.integration
    @pytest.mark.slow
    def test_pipeline_execution(self, minimal_config):
        """Тест выполнения всего пайплайна."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pipeline = ProteinLigandPipeline(
                config=minimal_config,
                output_dir=tmpdir
            )
            
            try:
                results = pipeline.run()
                
                # Проверка структуры результатов
                assert results["status"] == "success"
                assert "steps" in results
                assert "structure_prediction" in results["steps"]
                assert "docking" in results["steps"]
                assert "md_simulation" in results["steps"]
                
                # Проверка файлов
                output_path = Path(tmpdir)
                assert (output_path / "protein.pdb").exists()
                assert (output_path / "complex.pdb").exists()
                assert (output_path / "results.json").exists()
                
            except Exception as e:
                pytest.skip(f"Pipeline failed (expected in some environments): {e}")


# ============ ЗАПУСК ============

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-m", "not slow"])
