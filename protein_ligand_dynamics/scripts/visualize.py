#!/usr/bin/env python3
"""
Скрипт визуализации структур и траекторий.

Использование:
    python scripts/03_visualize.py --pdb results/complex.pdb \
                                    --trajectory results/trajectory.dcd \
                                    --output visualizations
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from utils.visualization import MolecularVisualizer
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def main():
    parser = argparse.ArgumentParser(description="Визуализация молекулярных структур")
    
    parser.add_argument("--pdb", help="Путь к PDB файлу")
    parser.add_argument("--trajectory", help="Путь к траектории")
    parser.add_argument("--output", default="visualizations", help="Папка для результатов")
    parser.add_argument("--create-gif", action="store_true", help="Создать GIF анимацию")
    parser.add_argument("--frame-skip", type=int, default=10, help="Пропуск кадров для GIF")
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    viz = MolecularVisualizer()
    
    print("="*70)
    print("  ВИЗУАЛИЗАЦИЯ")
    print("="*70)
    
    # === СТАТИЧЕСКИЕ ИЗОБРАЖЕНИЯ ===
    
    if args.pdb:
        print(f"\n[1] Обработка PDB: {args.pdb}")
        
        # Карта контактов
        print("    Генерация карты контактов...")
        viz.plot_contact_map(
            args.pdb,
            cutoff=8.0,
            output_file=str(output_dir / "contact_map.png")
        )
        
        # 3D структура (matplotlib)
        print("    3D представление...")
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        coords = []
        with open(args.pdb) as f:
            for line in f:
                if line.startswith("ATOM") and " CA " in line:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
        
        coords = np.array(coords)
        
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 
                'o-', markersize=3, linewidth=0.7, alpha=0.7)
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")
        ax.set_title("3D Структура белка (C-alpha trace)")
        
        plt.tight_layout()
        plt.savefig(output_dir / "structure_3d.png", dpi=300)
        print(f"    ✓ {output_dir / 'structure_3d.png'}")
    
    # === АНИМАЦИЯ ТРАЕКТОРИИ ===
    
    if args.trajectory and args.pdb:
        print(f"\n[2] Обработка траектории: {args.trajectory}")
        
        traj = md.load(args.trajectory, top=args.pdb)
        print(f"    Кадров: {traj.n_frames}")
        
        # Подвыборка
        traj_sub = traj[::args.frame_skip]
        print(f"    Подвыборка: {traj_sub.n_frames} кадров")
        
        if args.create_gif:
            print("    Создание GIF анимации...")
            
            # Выбор C-alpha
            ca_indices = traj_sub.topology.select("name CA")
            
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            
            def update(frame_idx):
                ax.clear()
                coords = traj_sub.xyz[frame_idx, ca_indices] * 10  # нм → Å
                
                ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
                       'o-', markersize=2, linewidth=0.5, alpha=0.7)
                
                ax.set_xlabel("X (Å)")
                ax.set_ylabel("Y (Å)")
                ax.set_zlabel("Z (Å)")
                ax.set_title(f"Кадр {frame_idx * args.frame_skip}")
                
                # Фиксированные пределы
                all_coords = traj_sub.xyz[:, ca_indices] * 10
                ax.set_xlim(all_coords[:, :, 0].min(), all_coords[:, :, 0].max())
                ax.set_ylim(all_coords[:, :, 1].min(), all_coords[:, :, 1].max())
                ax.set_zlim(all_coords[:, :, 2].min(), all_coords[:, :, 2].max())
            
            anim = FuncAnimation(fig, update, frames=traj_sub.n_frames, interval=100)
            
            gif_file = output_dir / "trajectory.gif"
            writer = PillowWriter(fps=10)
            anim.save(gif_file, writer=writer)
            
            print(f"    ✓ {gif_file}")
        
        # RMSD во времени с наложением структуры
        print("    RMSD эволюция с структурой...")
        
        rmsd = md.rmsd(traj, traj, 0) * 10  # Å
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # RMSD график
        ax1.plot(rmsd, linewidth=0.8)
        ax1.set_xlabel("Кадр")
        ax1.set_ylabel("RMSD (Å)")
        ax1.set_title("RMSD относительно начальной структуры")
        ax1.grid(alpha=0.3)
        
        # Наложение начальной и конечной структур
        ca_indices = traj.topology.select("name CA")
        
        coords_initial = traj.xyz[0, ca_indices] * 10
        coords_final = traj.xyz[-1, ca_indices] * 10
        
        ax2 = fig.add_subplot(122, projection='3d')
        ax2.plot(coords_initial[:, 0], coords_initial[:, 1], coords_initial[:, 2],
                'b-', linewidth=0.8, alpha=0.7, label='Начальная')
        ax2.plot(coords_final[:, 0], coords_final[:, 1], coords_final[:, 2],
                'r-', linewidth=0.8, alpha=0.7, label='Конечная')
        
        ax2.set_xlabel("X (Å)")
        ax2.set_ylabel("Y (Å)")
        ax2.set_zlabel("Z (Å)")
        ax2.set_title("Наложение структур")
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig(output_dir / "trajectory_overview.png", dpi=300)
        print(f"    ✓ {output_dir / 'trajectory_overview.png'}")
    
    # === ФИНАЛ ===
    
    print("\n" + "="*70)
    print("  ✓ ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА")
    print(f"  Результаты: {output_dir}/")
    print("="*70)


if __name__ == "__main__":
    main()
