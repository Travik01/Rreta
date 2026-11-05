import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from structure_prediction.esmfold_predictor import ESMFoldPredictor
from structure_prediction.quality_checker import classify_quality


class TestESMFoldPredictor:
    @pytest.fixture
    def predictor(self):
        return ESMFoldPredictor(use_local=False)

    def test_sequence_validation(self, predictor):
        assert predictor._validate_sequence("acd efg\nhi") == "ACDEFGHI"
        assert predictor._validate_sequence("ACDEFGHIKLMNPQRSTVWY") == "ACDEFGHIKLMNPQRSTVWY"

        with pytest.raises(ValueError):
            predictor._validate_sequence("ACDEFXYZ")
        with pytest.raises(ValueError):
            predictor._validate_sequence("")
        with pytest.raises(ValueError):
            predictor._validate_sequence("   \n  ")

    def test_plddt_extraction(self, predictor):
        def pdb_line(serial, name, x, y, z, bfactor):
            return (
                f"ATOM  {serial:5d} {name:>4s} MET A   1"
                f"   {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{bfactor:6.2f}           C"
            )

        pdb_string = "\n".join(
            [
                pdb_line(1, "N", 10.0, 10.0, 10.0, 95.50),
                pdb_line(2, "CA", 11.0, 11.0, 11.0, 92.30),
                pdb_line(3, "C", 12.0, 12.0, 12.0, 90.10),
            ]
        )
        scores = predictor._extract_plddt(pdb_string)
        assert len(scores) == 1
        assert scores[0] == 92.30

    @pytest.mark.slow
    def test_api_prediction(self, predictor):
        sequence = "MKFLKFSLLTAVLLSVVFAFSSCGDDDDTGYLPPS"
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = Path(tmpdir) / "test.pdb"
            pdb_string, metadata = predictor.predict(
                sequence,
                output_path=str(output_file),
            )

            assert output_file.exists()
            assert len(pdb_string) > 100
            assert metadata["method"] == "ESMFold_API"
            assert metadata["sequence_length"] == len(sequence)
            assert 0 <= metadata["plddt_mean"] <= 100


def test_plddt_quality_classification():
    assert classify_quality(95) == "Отлично"
    assert classify_quality(75) == "Хорошо"
    assert classify_quality(60) == "Средне"
    assert classify_quality(40) == "Плохо"
