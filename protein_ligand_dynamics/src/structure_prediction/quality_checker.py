"""
Проверка качества предсказанной структуры белка.

Вспомогательный модуль для ESMFoldPredictor:
классификация качества по pLDDT и базовые проверки структуры.
"""

from typing import Dict, List, Optional


def classify_quality(plddt: float) -> str:
    """
    Классификация качества структуры по среднему pLDDT.

    Args:
        plddt: Средний pLDDT (0-100)

    Returns:
        Строка с категорией качества
    """
    if plddt > 90:
        return "Отлично"
    elif plddt > 70:
        return "Хорошо"
    elif plddt > 50:
        return "Средне"
    else:
        return "Плохо"


class QualityChecker:
    """
    Проверка качества предсказанной структуры.

    Основной метод: classify(plddt_scores) -> dict с оценками.
    """

    # Пороги pLDDT (по рекомендациям AlphaFold/ESMFold)
    THRESHOLDS = {
        "very_high": 90,  # очень высокая уверенность
        "high": 70,       # хорошая модель
        "low": 50,        # низкая уверенность
    }

    def __init__(self, plddt_scores: Optional[List[float]] = None):
        self.plddt_scores = plddt_scores or []

    @property
    def mean_plddt(self) -> float:
        """Средний pLDDT."""
        if not self.plddt_scores:
            return 0.0
        return sum(self.plddt_scores) / len(self.plddt_scores)

    def fraction_above(self, threshold: float) -> float:
        """Доля остатков с pLDDT выше порога (0-1)."""
        if not self.plddt_scores:
            return 0.0
        return sum(1 for s in self.plddt_scores if s > threshold) / len(self.plddt_scores)

    def classify(self) -> Dict:
        """
        Полная оценка качества.

        Returns:
            dict с полями: mean_plddt, quality, fraction_very_high, ...
        """
        mean = self.mean_plddt

        result = {
            "mean_plddt": mean,
            "quality": classify_quality(mean),
            "n_residues": len(self.plddt_scores),
            "fraction_very_high": self.fraction_above(self.THRESHOLDS["very_high"]),
            "fraction_high": self.fraction_above(self.THRESHOLDS["high"]),
            "fraction_low": self.fraction_above(self.THRESHOLDS["low"]),
        }

        return result


# ============ ТЕСТ ============
if __name__ == "__main__":
    checker = QualityChecker([95.0, 88.0, 72.0, 45.0, 91.0])
    for key, value in checker.classify().items():
        print(f"  {key}: {value}")
