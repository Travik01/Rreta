import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from full_pipeline import ProteinLigandPipeline


@pytest.mark.integration
@pytest.mark.slow
def test_pipeline_execution():
    config = {
        "structure_prediction": {
            "sequence": "MKFLKFSLLTAVLLSVVFAFSSCG",
            "use_local": False,
        },
        "docking": {
            "ligand_smiles": "CC(=O)O",
            "num_poses": 2,
        },
        "md_simulation": {
            "use_mace": False,
            "use_gpu": False,
            "temperature": 300.0,
            "num_steps": 100,
            "report_interval": 50,
        },
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        results = ProteinLigandPipeline(config, tmpdir).run()
        output_path = Path(tmpdir)

        assert results["status"] == "success"
        assert "structure_prediction" in results["steps"]
        assert "docking" in results["steps"]
        assert "md_simulation" in results["steps"]
        assert (output_path / "protein.pdb").exists()
        assert (output_path / "complex.pdb").exists()
        assert (output_path / "results.json").exists()
