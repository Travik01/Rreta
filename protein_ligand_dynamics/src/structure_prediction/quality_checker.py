"""Оценка качества структуры по pLDDT."""

from typing import ClassVar


def classify_quality(plddt: float) -> str:
    if plddt > 90:
        return "Отлично"
    if plddt > 70:
        return "Хорошо"
    if plddt > 50:
        return "Средне"
    return "Плохо"


class QualityChecker:
    THRESHOLDS: ClassVar[dict[str, int]] = {"very_high": 90, "high": 70, "low": 50}

    def __init__(self, plddt_scores: list[float] | None = None):
        self.plddt_scores = plddt_scores or []

    @property
    def mean_plddt(self) -> float:
        if not self.plddt_scores:
            return 0.0
        return sum(self.plddt_scores) / len(self.plddt_scores)

    def fraction_above(self, threshold: float) -> float:
        if not self.plddt_scores:
            return 0.0
        return sum(score > threshold for score in self.plddt_scores) / len(self.plddt_scores)

    def classify(self) -> dict:
        return {
            "mean_plddt": self.mean_plddt,
            "quality": classify_quality(self.mean_plddt),
            "n_residues": len(self.plddt_scores),
            "fraction_very_high": self.fraction_above(self.THRESHOLDS["very_high"]),
            "fraction_high": self.fraction_above(self.THRESHOLDS["high"]),
            "fraction_low": self.fraction_above(self.THRESHOLDS["low"]),
        }
