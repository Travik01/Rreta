#!/usr/bin/env python3
"""
Скрипт для анализа результатов MD-симуляции.

Использование:
    python scripts/02_analyze_results.py --trajectory results/trajectory.dcd \
                                         --topology results/complex.pdb \
                                         --output analysis_report
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from analysis.trajectory_analyzer import TrajectoryAnalyzer
from utils.visualization import MolecularVisualizer
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Анализ MD-траектории")
    
    parser.add_argument(
        "--trajectory",
        required=True,
        help="Путь к траектории (.dcd, .xtc, .traj)"
    )
    
    parser.add_argument(
        "--topology",
        required=True,
        help="Путь к топологии (.pdb)"
    )
    
    parser.add_argument(
        "--output",
        default="analysis_output",
        help="Папка для результатов"
    )
    
    parser.add_argument(
        "--rmsd-selection",
        default="protein and name CA",
        help="Селекция для RMSD"
    )
    
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=4.5,
        help="Cutoff для контактов (Å)"
    )
    
    args = parser.parse_args()
    
    print("="*70)
    print("  АНАЛИЗ MD-ТРАЕКТОРИИ")
    print("="*70)
    
    # Загрузка данных
    analyzer = TrajectoryAnalyzer(
        trajectory_file=args.trajectory,
        topology_file=args.topology
    )
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # === БАЗОВЫЕ МЕТРИКИ ===
    
    print("\n[1] RMSD...")
    rmsd = analyzer.compute_rmsd(selection=args.rmsd_selection)
    
    print("\n[2] RMSF...")
    rmsf, residues = analyzer.compute_rmsf()
    
    print("\n[3] Radius of Gyration...")
    rg = analyzer.compute_radius_of_gyration()
    
    print("\n[4] Контакты...")
    contact_freq, contact_pairs = analyzer.compute_contacts(
        cutoff=args.contact_cutoff / 10  # Å → нм
    )
    
    print("\n[5] Расстояние белок-лиганд...")
    try:
        distances = analyzer.compute_ligand_protein_distance()
    except ValueError as e:
        print(f"    Пропуск: {e}")
        distances = None
    
    # === РАСШИРЕННЫЙ АНАЛИЗ ===
    
    print("\n[6] Дополнительная статистика...")
    
    # Конвергенция RMSD
    rmsd_first_half = rmsd[:len(rmsd)//2]
    rmsd_second_half = rmsd[len(rmsd)//2:]
    
    convergence = {
        "rmsd_first_half_mean": np.mean(rmsd_first_half),
        "rmsd_second_half_mean": np.mean(rmsd_second_half),
        "rmsd_drift": np.mean(rmsd_second_half) - np.mean(rmsd_first_half)
    }
    
    print(f"    RMSD drift: {convergence['rmsd_drift']:.2f} Å")
    
    if abs(convergence['rmsd_drift']) < 0.5:
        print("    ✓ Система стабильна (drift < 0.5 Å)")
    else:
        print("    ⚠ Система может быть нестабильной")
    
    # === ВИЗУАЛИЗАЦИЯ ===
    
    print("\n[7] Генерация графиков...")
    
    # Комбинированный график
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # RMSD
    axes[0, 0].plot(rmsd, linewidth=0.8)
    axes[0, 0].set_xlabel("Кадр")
    axes[0, 0].set_ylabel("RMSD (Å)")
    axes[0, 0].set_title("RMSD во времени")
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].axhline(np.mean(rmsd), color='r', linestyle='--', alpha=0.5, label='Среднее')
    axes[0, 0].legend()
    
    # RMSF
    axes[0, 1].bar(residues, rmsf, width=1.0, alpha=0.7)
    axes[0, 1].set_xlabel("Номер остатка")
    axes[0, 1].set_ylabel("RMSF (Å)")
    axes[0, 1].set_title("Флуктуации остатков")
    axes[0, 1].grid(alpha=0.3)
    
    # Rg
    axes[1, 0].plot(rg, linewidth=0.8, color='green')
    axes[1, 0].set_xlabel("Кадр")
    axes[1, 0].set_ylabel("Rg (Å)")
    axes[1, 0].set_title("Радиус гирации")
    axes[1, 0].grid(alpha=0.3)
    axes[1, 0].axhline(np.mean(rg), color='r', linestyle='--', alpha=0.5)
    
    # Расстояние
    if distances is not None:
        axes[1, 1].plot(distances, linewidth=0.8, color='purple')
        axes[1, 1].set_xlabel("Кадр")
        axes[1, 1].set_ylabel("Расстояние (Å)")
        axes[1, 1].set_title("Белок-лиганд расстояние")
        axes[1, 1].grid(alpha=0.3)
        axes[1, 1].axhline(np.mean(distances), color='r', linestyle='--', alpha=0.5)
    else:
        axes[1, 1].text(0.5, 0.5, "Лиганд не найден", 
                       ha='center', va='center', transform=axes[1, 1].transAxes)
        axes[1, 1].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_dir / "combined_analysis.png", dpi=300, bbox_inches='tight')
    print(f"    ✓ {output_dir / 'combined_analysis.png'}")
    
    # Гистограмма RMSD
    plt.figure(figsize=(8, 5))
    plt.hist(rmsd, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel("RMSD (Å)")
    plt.ylabel("Частота")
    plt.title("Распределение RMSD")
    plt.axvline(np.mean(rmsd), color='r', linestyle='--', label=f'Среднее: {np.mean(rmsd):.2f} Å')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "rmsd_histogram.png", dpi=300)
    print(f"    ✓ {output_dir / 'rmsd_histogram.png'}")
    
    # === СОХРАНЕНИЕ ДАННЫХ ===
    
    print("\n[8] Сохранение данных...")
    
    import pandas as pd
    
    df = pd.DataFrame({
        "frame": range(len(rmsd)),
        "rmsd": rmsd,
        "rg": rg,
    })
    
    if distances is not None:
        df["ligand_distance"] = distances
    
    csv_file = output_dir / "metrics.csv"
    df.to_csv(csv_file, index=False)
    print(f"    ✓ {csv_file}")
    
    # Сводка в JSON
    import json
    
    summary = {
        "num_frames": len(rmsd),
        "rmsd": {
            "mean": float(np.mean(rmsd)),
            "std": float(np.std(rmsd)),
            "min": float(np.min(rmsd)),
            "max": float(np.max(rmsd))
        },
        "rmsf": {
            "mean": float(np.mean(rmsf)),
            "max": float(np.max(rmsf)),
            "max_residue": int(residues[np.argmax(rmsf)])
        },
        "rg": {
            "mean": float(np.mean(rg)),
            "std": float(np.std(rg))
        },
        "convergence": convergence
    }
    
    if distances is not None:
        summary["ligand_distance"] = {
            "mean": float(np.mean(distances)),
            "min": float(np.min(distances))
        }
    
    json_file = output_dir / "summary.json"
    with open(json_file, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"    ✓ {json_file}")
    
    # === ФИНАЛ ===
    
    print("\n" + "="*70)
    print("  ✓ АНАЛИЗ ЗАВЕРШЁН")
    print(f"  Результаты: {output_dir}/")
    print("="*70)
    print("\nКлючевые метрики:")
    print(f"  • RMSD: {summary['rmsd']['mean']:.2f} ± {summary['rmsd']['std']:.2f} Å")
    print(f"  • Rg: {summary['rg']['mean']:.2f} ± {summary['rg']['std']:.2f} Å")
    
    if distances is not None:
        print(f"  • Мин. расстояние до лиганда: {summary['ligand_distance']['min']:.2f} Å")
    
    print(f"  • Самый подвижный остаток: {summary['rmsf']['max_residue']} (RMSF={summary['rmsf']['max']:.2f} Å)")


if __name__ == "__main__":
    main()
