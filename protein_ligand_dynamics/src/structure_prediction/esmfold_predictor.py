"""Получение PDB-структуры белка через ESMFold."""

import logging
import warnings
from pathlib import Path

import requests

from .quality_checker import classify_quality

try:
    import torch
except ImportError:
    torch = None

logger = logging.getLogger(__name__)


class ESMFoldPredictor:
    API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    def __init__(
        self,
        use_local: bool = False,
        device: str = "auto",
        cache_dir: str = "./models",
    ):
        self.use_local = use_local
        self.model = None
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        if use_local and torch is None:
            raise ImportError("Для локального ESMFold установите PyTorch: pip install torch")

        if device == "auto":
            self.device = "cuda" if (torch is not None and torch.cuda.is_available()) else "cpu"
        else:
            self.device = device

        if use_local:
            self._load_local_model()

    def _load_local_model(self):
        logger.info("Загрузка локальной модели ESMFold")
        try:
            import esm
        except ImportError as exc:
            raise ImportError("Для локального ESMFold установите: pip install fair-esm") from exc

        self.model = esm.pretrained.esmfold_v1().eval()
        if self.device == "cuda":
            self.model = self.model.cuda()
            logger.info("Локальная модель на GPU: %s", torch.cuda.get_device_name(0))
        else:
            logger.info("Локальная модель на CPU")
        self.model.set_chunk_size(64)

    def predict(
        self,
        sequence: str,
        output_path: str | None = None,
        return_confidence: bool = True,
    ) -> tuple[str, dict]:
        sequence = self._validate_sequence(sequence)
        logger.info("Предсказание структуры, %d а.о.", len(sequence))

        if self.use_local:
            pdb_string, metadata = self._predict_local(sequence)
        else:
            pdb_string, metadata = self._predict_api(sequence)

        if output_path:
            path = Path(output_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(pdb_string, encoding="utf-8")
            metadata["output_file"] = str(path)
            logger.info("Структура сохранена: %s", path)

        if return_confidence:
            scores = self._extract_plddt(pdb_string)
            metadata["plddt_scores"] = scores
            metadata["plddt_mean"] = sum(scores) / len(scores) if scores else 0.0
            metadata["quality"] = classify_quality(metadata["plddt_mean"])
            logger.info(
                "Средний pLDDT: %.1f (%s)",
                metadata["plddt_mean"],
                metadata["quality"],
            )

        return pdb_string, metadata

    def _predict_api(self, sequence: str) -> tuple[str, dict]:
        try:
            response = requests.post(
                self.API_URL,
                data=sequence,
                headers={"Content-Type": "text/plain"},
                timeout=300,
            )
        except requests.Timeout as exc:
            raise RuntimeError("ESMFold API не ответил за 300 с") from exc

        if response.status_code == 200:
            return response.text, {
                "method": "ESMFold_API",
                "sequence_length": len(sequence),
                "status": "success",
            }
        if response.status_code == 429:
            raise RuntimeError(
                "ESMFold API отклонил запрос (429); повторите позже или используйте локальную модель"
            )
        raise RuntimeError(f"ESMFold API вернул {response.status_code}: {response.text}")

    def _predict_local(self, sequence: str) -> tuple[str, dict]:
        with torch.no_grad():
            pdb_string = self.model.infer_pdb(sequence)

        return pdb_string, {
            "method": "ESMFold_local",
            "sequence_length": len(sequence),
            "device": self.device,
            "status": "success",
        }

    @staticmethod
    def _validate_sequence(sequence: str) -> str:
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        sequence = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")

        if not sequence:
            raise ValueError("Пустая последовательность")

        invalid = sorted(set(sequence) - valid_aa)
        if invalid:
            raise ValueError(f"Недопустимые символы: {''.join(invalid)}")

        if len(sequence) > 1000:
            warnings.warn(
                f"Длинная последовательность ({len(sequence)} а.о.): "
                "предсказание может быть медленным или неточным",
                stacklevel=2,
            )

        return sequence

    @staticmethod
    def _extract_plddt(pdb_string: str) -> list:
        scores = []
        for line in pdb_string.split("\n"):
            if line.startswith("ATOM") and len(line) >= 66:
                if line[12:16].strip() != "CA":
                    continue
                try:
                    scores.append(float(line[60:66].strip()))
                except ValueError:
                    continue
        return scores
