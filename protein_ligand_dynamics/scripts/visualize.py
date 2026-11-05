#!/usr/bin/env python3
"""Визуализация структур и траекторий."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from utils.visualization import MolecularVisualizer


def main():
    parser = argparse.ArgumentParser(description="Визуализация молекулярных структур")
    parser.add_argument("--pdb", help="Путь к PDB-файлу")
    parser.add_argument("--trajectory", help="Путь к траектории")
    parser.add_argument("--output", default="visualizations", help="Каталог результатов")
    parser.add_argument("--create-gif", action="store_true", help="Создать GIF")
    parser.add_argument(
        "--frame-skip",
        type=int,
        default=10,
        help="Шаг по кадрам для GIF",
    )
    args = parser.parse_args()

    if args.frame_skip < 1:
        parser.error("--frame-skip должен быть не меньше 1")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    viz = MolecularVisualizer()

    if args.pdb:
        contact_map = output_dir / "contact_map.png"
        viz.plot_contact_map(args.pdb, cutoff=8.0, output_file=str(contact_map))
        print(f"Сохранено: {contact_map}")

        coords = []
        with open(args.pdb, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("ATOM") and " CA " in line:
                    coords.append(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    )

        coords = np.asarray(coords)
        if coords.size:
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection="3d")
            ax.plot(
                coords[:, 0],
                coords[:, 1],
                coords[:, 2],
                "o-",
                markersize=3,
                linewidth=0.7,
                alpha=0.7,
            )
            ax.set_xlabel("X (Å)")
            ax.set_ylabel("Y (Å)")
            ax.set_zlabel("Z (Å)")
            ax.set_title("C-alpha trace")
            fig.tight_layout()
            structure_3d = output_dir / "structure_3d.png"
            fig.savefig(structure_3d, dpi=300)
            plt.close(fig)
            print(f"Сохранено: {structure_3d}")

    if args.trajectory and args.pdb:
        traj = md.load(args.trajectory, top=args.pdb)
        print(f"Кадров: {traj.n_frames}")
        traj_sub = traj[:: args.frame_skip]

        if args.create_gif:
            ca_indices = traj_sub.topology.select("name CA")
            all_coords = traj_sub.xyz[:, ca_indices] * 10

            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection="3d")

            def update(frame_idx):
                ax.clear()
                frame_coords = all_coords[frame_idx]
                ax.plot(
                    frame_coords[:, 0],
                    frame_coords[:, 1],
                    frame_coords[:, 2],
                    "o-",
                    markersize=2,
                    linewidth=0.5,
                    alpha=0.7,
                )
                ax.set_xlabel("X (Å)")
                ax.set_ylabel("Y (Å)")
                ax.set_zlabel("Z (Å)")
                ax.set_title(f"Кадр {frame_idx * args.frame_skip}")
                ax.set_xlim(all_coords[:, :, 0].min(), all_coords[:, :, 0].max())
                ax.set_ylim(all_coords[:, :, 1].min(), all_coords[:, :, 1].max())
                ax.set_zlim(all_coords[:, :, 2].min(), all_coords[:, :, 2].max())

            animation = FuncAnimation(
                fig, update, frames=traj_sub.n_frames, interval=100
            )
            gif_file = output_dir / "trajectory.gif"
            animation.save(gif_file, writer=PillowWriter(fps=10))
            plt.close(fig)
            print(f"Сохранено: {gif_file}")

        rmsd = md.rmsd(traj, traj, 0) * 10
        ca_indices = traj.topology.select("name CA")
        initial = traj.xyz[0, ca_indices] * 10
        final = traj.xyz[-1, ca_indices] * 10

        fig = plt.figure(figsize=(14, 5))
        ax1 = fig.add_subplot(121)
        ax1.plot(rmsd, linewidth=0.8)
        ax1.set_xlabel("Кадр")
        ax1.set_ylabel("RMSD (Å)")
        ax1.set_title("RMSD относительно первого кадра")
        ax1.grid(alpha=0.3)

        ax2 = fig.add_subplot(122, projection="3d")
        ax2.plot(
            initial[:, 0],
            initial[:, 1],
            initial[:, 2],
            "b-",
            linewidth=0.8,
            alpha=0.7,
            label="Начальная",
        )
        ax2.plot(
            final[:, 0],
            final[:, 1],
            final[:, 2],
            "r-",
            linewidth=0.8,
            alpha=0.7,
            label="Конечная",
        )
        ax2.set_xlabel("X (Å)")
        ax2.set_ylabel("Y (Å)")
        ax2.set_zlabel("Z (Å)")
        ax2.set_title("Наложение структур")
        ax2.legend()

        fig.tight_layout()
        overview = output_dir / "trajectory_overview.png"
        fig.savefig(overview, dpi=300)
        plt.close(fig)
        print(f"Сохранено: {overview}")


if __name__ == "__main__":
    main()
