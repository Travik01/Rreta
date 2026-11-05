#!/usr/bin/env python3
"""Анализ траектории MD."""

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from analysis.trajectory_analyzer import TrajectoryAnalyzer


def main():
    parser = argparse.ArgumentParser(description="Анализ MD-траектории")
    parser.add_argument("--trajectory", required=True, help="Путь к траектории (.dcd, .xtc)")
    parser.add_argument("--topology", required=True, help="Путь к топологии (.pdb)")
    parser.add_argument("--output", default="analysis_output", help="Каталог результатов")
    parser.add_argument(
        "--rmsd-selection",
        default="protein and name CA",
        help="Селекция для RMSD",
    )
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=4.5,
        help="Порог контактов, Å",
    )
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    analyzer = TrajectoryAnalyzer(args.trajectory, args.topology)

    rmsd = analyzer.compute_rmsd(selection=args.rmsd_selection)
    rmsf, residues = analyzer.compute_rmsf()
    rg = analyzer.compute_radius_of_gyration()
    contact_freq, contact_pairs = analyzer.compute_contacts(
        cutoff=args.contact_cutoff / 10
    )

    try:
        distances = analyzer.compute_ligand_protein_distance()
    except ValueError as exc:
        print(f"Расстояние до лиганда не рассчитано: {exc}")
        distances = None

    half = len(rmsd) // 2
    convergence = {
        "rmsd_first_half_mean": float(np.mean(rmsd[:half])),
        "rmsd_second_half_mean": float(np.mean(rmsd[half:])),
        "rmsd_drift": float(np.mean(rmsd[half:]) - np.mean(rmsd[:half])),
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    axes[0, 0].plot(rmsd, linewidth=0.8)
    axes[0, 0].set_xlabel("Кадр")
    axes[0, 0].set_ylabel("RMSD (Å)")
    axes[0, 0].set_title("RMSD")
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].axhline(
        np.mean(rmsd), color="r", linestyle="--", alpha=0.5, label="Среднее"
    )
    axes[0, 0].legend()

    axes[0, 1].bar(residues, rmsf, width=1.0, alpha=0.7)
    axes[0, 1].set_xlabel("Номер остатка")
    axes[0, 1].set_ylabel("RMSF (Å)")
    axes[0, 1].set_title("Флуктуации остатков")
    axes[0, 1].grid(alpha=0.3)

    axes[1, 0].plot(rg, linewidth=0.8, color="green")
    axes[1, 0].set_xlabel("Кадр")
    axes[1, 0].set_ylabel("Rg (Å)")
    axes[1, 0].set_title("Радиус гирации")
    axes[1, 0].grid(alpha=0.3)
    axes[1, 0].axhline(np.mean(rg), color="r", linestyle="--", alpha=0.5)

    if distances is not None:
        axes[1, 1].plot(distances, linewidth=0.8, color="purple")
        axes[1, 1].set_xlabel("Кадр")
        axes[1, 1].set_ylabel("Расстояние (Å)")
        axes[1, 1].set_title("Белок–лиганд")
        axes[1, 1].grid(alpha=0.3)
        axes[1, 1].axhline(
            np.mean(distances), color="r", linestyle="--", alpha=0.5
        )
    else:
        axes[1, 1].text(
            0.5,
            0.5,
            "Лиганд не найден",
            ha="center",
            va="center",
            transform=axes[1, 1].transAxes,
        )
        axes[1, 1].axis("off")

    fig.tight_layout()
    combined = output_dir / "combined_analysis.png"
    fig.savefig(combined, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Сохранено: {combined}")

    plt.figure(figsize=(8, 5))
    plt.hist(rmsd, bins=50, alpha=0.7, edgecolor="black")
    plt.xlabel("RMSD (Å)")
    plt.ylabel("Частота")
    plt.title("Распределение RMSD")
    plt.axvline(
        np.mean(rmsd),
        color="r",
        linestyle="--",
        label=f"Среднее: {np.mean(rmsd):.2f} Å",
    )
    plt.legend()
    plt.tight_layout()
    histogram = output_dir / "rmsd_histogram.png"
    plt.savefig(histogram, dpi=300)
    plt.close()
    print(f"Сохранено: {histogram}")

    df = pd.DataFrame({"frame": range(len(rmsd)), "rmsd": rmsd, "rg": rg})
    if distances is not None:
        df["ligand_distance"] = distances
    csv_file = output_dir / "metrics.csv"
    df.to_csv(csv_file, index=False)
    print(f"Сохранено: {csv_file}")

    summary = {
        "num_frames": len(rmsd),
        "rmsd": {
            "mean": float(np.mean(rmsd)),
            "std": float(np.std(rmsd)),
            "min": float(np.min(rmsd)),
            "max": float(np.max(rmsd)),
        },
        "rmsf": {
            "mean": float(np.mean(rmsf)),
            "max": float(np.max(rmsf)),
            "max_residue": int(residues[np.argmax(rmsf)]),
        },
        "rg": {
            "mean": float(np.mean(rg)),
            "std": float(np.std(rg)),
        },
        "contacts": {
            "num_residues": len(contact_pairs),
            "mean_frequency": float(np.mean(contact_freq)) if contact_freq.size else 0.0,
        },
        "convergence": convergence,
    }

    if distances is not None:
        summary["ligand_distance"] = {
            "mean": float(np.mean(distances)),
            "min": float(np.min(distances)),
        }

    json_file = output_dir / "summary.json"
    with open(json_file, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
    print(f"Сохранено: {json_file}")

    print(f"RMSD: {summary['rmsd']['mean']:.2f} ± {summary['rmsd']['std']:.2f} Å")
    print(f"Rg: {summary['rg']['mean']:.2f} ± {summary['rg']['std']:.2f} Å")
    if distances is not None:
        print(
            f"Минимальное расстояние до лиганда: "
            f"{summary['ligand_distance']['min']:.2f} Å"
        )
    print(
        f"Самый подвижный остаток: {summary['rmsf']['max_residue']} "
        f"(RMSF={summary['rmsf']['max']:.2f} Å)"
    )
    if abs(convergence["rmsd_drift"]) < 0.5:
        print("Дрейф RMSD < 0.5 Å")
    else:
        print("Дрейф RMSD >= 0.5 Å")


if __name__ == "__main__":
    main()
