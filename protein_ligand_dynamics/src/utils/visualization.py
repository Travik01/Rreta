"""Визуализация структур и траекторий."""

import io
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from scipy.spatial.distance import cdist

try:
    import py3Dmol
except ImportError:
    py3Dmol = None

logger = logging.getLogger(__name__)


class MolecularVisualizer:
    @staticmethod
    def show_structure(
        pdb_file: str,
        style: str = "cartoon",
        color_scheme: str = "spectrum",
        width: int = 800,
        height: int = 600,
        surface_opacity: float = 0.7,
    ):
        if py3Dmol is None:
            raise ImportError("Для 3D-визуализации установите: pip install py3Dmol")

        pdb_string = Path(pdb_file).read_text(encoding="utf-8")
        view = py3Dmol.view(width=width, height=height)
        view.addModel(pdb_string, "pdb")
        view.setStyle({"chain": "A"}, {style: {"color": color_scheme}})
        view.setStyle(
            {"chain": "L"},
            {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}},
        )
        view.addSurface(
            py3Dmol.VDW,
            {"opacity": surface_opacity, "color": "white"},
            {"chain": "A"},
        )
        view.zoomTo()
        return view

    @staticmethod
    def show_trajectory(
        trajectory_file: str,
        topology_file: str,
        frame_skip: int = 10,
        width: int = 800,
        height: int = 600,
    ):
        if py3Dmol is None:
            raise ImportError("Для 3D-визуализации установите: pip install py3Dmol")

        traj = md.load(trajectory_file, top=topology_file)[::frame_skip]
        view = py3Dmol.view(width=width, height=height)

        for frame in traj:
            buffer = io.StringIO()
            frame.save_pdb(buffer)
            view.addModel(buffer.getvalue(), "pdb")

        view.setStyle({"cartoon": {"color": "spectrum"}})
        view.animate({"loop": "forward", "interval": 100})
        view.zoomTo()
        return view

    @staticmethod
    def create_comparison_plot(
        original_pdb: str,
        simulated_pdb: str,
        output_file: str = "comparison.png",
    ):
        fig = plt.figure(figsize=(15, 7))
        ax1 = fig.add_subplot(121, projection="3d")
        MolecularVisualizer._plot_structure_3d(original_pdb, ax1, "Начальная структура")
        ax2 = fig.add_subplot(122, projection="3d")
        MolecularVisualizer._plot_structure_3d(simulated_pdb, ax2, "После симуляции")
        fig.tight_layout()
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        logger.info("Сравнение сохранено: %s", output_file)

    @staticmethod
    def _plot_structure_3d(pdb_file: str, ax, title: str):
        coords = []
        with open(pdb_file, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("ATOM") and " CA " in line:
                    coords.append(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                    )

        coords = np.asarray(coords)
        if coords.size:
            ax.plot(
                coords[:, 0],
                coords[:, 1],
                coords[:, 2],
                "o-",
                markersize=2,
                linewidth=0.5,
                alpha=0.6,
            )

        ax.set_title(title)
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")

    @staticmethod
    def plot_contact_map(
        pdb_file: str,
        cutoff: float = 8.0,
        output_file: str = "contact_map.png",
    ):
        traj = md.load(pdb_file)
        ca_indices = traj.topology.select("name CA")
        distances = cdist(traj.xyz[0, ca_indices], traj.xyz[0, ca_indices])
        contact_matrix = (distances < cutoff / 10).astype(int)

        plt.figure(figsize=(10, 8))
        plt.imshow(contact_matrix, cmap="Blues", interpolation="nearest")
        plt.colorbar(label="Контакт")
        plt.xlabel("Остаток")
        plt.ylabel("Остаток")
        plt.title(f"Карта контактов (cutoff={cutoff} Å)")
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        logger.info("Карта контактов сохранена: %s", output_file)
