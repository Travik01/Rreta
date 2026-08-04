"""
ESMFold predictor с API и локальной версией.
"""

import requests
import warnings
from pathlib import Path
from typing import Optional, Dict, Tuple

# torch нужен ТОЛЬКО для локальной версии ESMFold.
# Для API-режима он не требуется, поэтому импортируем его лениво,
# чтобы модуль работал без установки тяжёлого PyTorch.
try:
    import torch
    _TORCH_AVAILABLE = True
except ImportError:
    torch = None
    _TORCH_AVAILABLE = False


class ESMFoldPredictor:
    """
    Предсказание 3D структуры белка через ESMFold.

    Методы:
    - API (бесплатно, без ключа, до 800 а.о.)
    - Локально (требует GPU ~16GB)
    """

    def __init__(
        self,
        use_local: bool = False,
        device: str = "auto",
        cache_dir: str = "./models"
    ):
        self.use_local = use_local
        self.model = None
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        if use_local and not _TORCH_AVAILABLE:
            raise ImportError(
                "Для локального ESMFold установите PyTorch: pip install torch"
            )

        if device == "auto":
            self.device = "cuda" if (_TORCH_AVAILABLE and torch.cuda.is_available()) else "cpu"
        else:
            self.device = device

        if use_local:
            self._load_local_model()

    def _load_local_model(self):
        """Загрузка локального ESMFold."""
        print("[ESMFold] Загрузка локальной модели...")

        try:
            import esm

            # Загрузка модели
            self.model = esm.pretrained.esmfold_v1()
            self.model = self.model.eval()

            if self.device == "cuda":
                self.model = self.model.cuda()
                print(f"           GPU: {torch.cuda.get_device_name(0)}")
            else:
                print("           CPU mode (медленно!)")

            # Оптимизация памяти
            self.model.set_chunk_size(64)

            print("[ESMFold] ✓ Модель загружена")

        except ImportError:
            raise ImportError(
                "Для локального ESMFold установи: pip install fair-esm"
            )

    def predict(
        self,
        sequence: str,
        output_path: Optional[str] = None,
        return_confidence: bool = True
    ) -> Tuple[str, Dict]:
        """
        Предсказывает структуру.

        Args:
            sequence: Аминокислотная последовательность
            output_path: Путь для сохранения PDB
            return_confidence: Вернуть pLDDT scores

        Returns:
            (pdb_string, metadata)
        """
        # Валидация и нормализация последовательности
        # (убираем пробелы/переводы строк, приводим к верхнему регистру)
        sequence = self._validate_sequence(sequence)

        print(f"[ESMFold] Предсказание структуры ({len(sequence)} а.о.)...")

        if self.use_local:
            pdb_string, metadata = self._predict_local(sequence)
        else:
            pdb_string, metadata = self._predict_api(sequence)

        # Сохранение
        if output_path:
            output_path = Path(output_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with open(output_path, "w") as f:
                f.write(pdb_string)

            metadata["output_file"] = str(output_path)
            print(f"[ESMFold] ✓ Сохранено: {output_path}")

        # Извлечение pLDDT
        if return_confidence:
            plddt_scores = self._extract_plddt(pdb_string)
            metadata["plddt_mean"] = sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0
            metadata["plddt_scores"] = plddt_scores

            print(f"[ESMFold] ✓ Средний pLDDT: {metadata['plddt_mean']:.1f}")

            # Интерпретация
            if metadata["plddt_mean"] > 90:
                quality = "Отлично"
            elif metadata["plddt_mean"] > 70:
                quality = "Хорошо"
            elif metadata["plddt_mean"] > 50:
                quality = "Средне (осторожно)"
            else:
                quality = "Плохо (не используй)"

            print(f"           Качество: {quality}")

        return pdb_string, metadata

    def _predict_api(self, sequence: str) -> Tuple[str, Dict]:
        """Предсказание через ESMAtlas API."""
        api_url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

        try:
            response = requests.post(
                api_url,
                data=sequence,
                headers={"Content-Type": "text/plain"},
                timeout=300  # 5 минут
            )

            if response.status_code == 200:
                pdb_string = response.text
                metadata = {
                    "method": "ESMFold_API",
                    "sequence_length": len(sequence),
                    "status": "success"
                }
                return pdb_string, metadata

            elif response.status_code == 429:
                raise RuntimeError(
                    "API rate limit. Подожди 1 минуту или используй локальную версию."
                )

            else:
                raise RuntimeError(
                    f"API error {response.status_code}: {response.text}"
                )

        except requests.Timeout:
            raise RuntimeError(
                "API timeout. Попробуй снова или используй локальную версию."
            )

    def _predict_local(self, sequence: str) -> Tuple[str, Dict]:
        """Локальное предсказание."""
        with torch.no_grad():
            pdb_string = self.model.infer_pdb(sequence)

        metadata = {
            "method": "ESMFold_local",
            "sequence_length": len(sequence),
            "device": self.device,
            "status": "success"
        }

        return pdb_string, metadata

    @staticmethod
    def _validate_sequence(sequence: str) -> str:
        """
        Проверка и нормализация валидности последовательности.

        Возвращает очищенную последовательность (верхний регистр,
        без пробелов и переводов строк).
        """
        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        sequence = sequence.upper().replace(" ", "").replace("\n", "").replace("\r", "")

        if not sequence:
            raise ValueError("Пустая последовательность")

        invalid_chars = set(sequence) - valid_aa
        if invalid_chars:
            raise ValueError(
                f"Невалидные символы: {invalid_chars}\n"
                f"Допустимые: {valid_aa}"
            )

        if len(sequence) > 1000:
            warnings.warn(
                f"Длинная последовательность ({len(sequence)} а.о.). "
                f"Предсказание может быть медленным или неточным."
            )

        return sequence

    @staticmethod
    def _extract_plddt(pdb_string: str) -> list:
        """Извлекает pLDDT confidence scores из B-factor (колонки 61-66)."""
        scores = []
        for line in pdb_string.split("\n"):
            if line.startswith("ATOM") and len(line) >= 66:
                atom_name = line[12:16].strip()
                if atom_name == "CA":  # Только C-alpha
                    try:
                        b_factor = float(line[60:66].strip())
                        scores.append(b_factor)
                    except (ValueError, IndexError):
                        pass
        return scores


# ============ ТЕСТ ============
if __name__ == "__main__":
    # Тест с API
    predictor = ESMFoldPredictor(use_local=False)

    # Короткий тестовый белок (WW domain)
    test_seq = "GSKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEP"

    pdb_data, metadata = predictor.predict(
        test_seq,
        output_path="test_output/protein.pdb"
    )

    print("\nМетаданные:")
    for k, v in metadata.items():
        if k != "plddt_scores":
            print(f"  {k}: {v}")
