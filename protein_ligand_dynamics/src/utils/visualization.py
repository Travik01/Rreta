"""
3D визуализация структур и траекторий.
"""

import io
import numpy as np
from pathlib import Path
from typing import Optional, List
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mdtraj as md
from scipy.spatial.distance import cdist

# py3Dmol загружается лениво (нужен только в Jupyter)
try:
    import py3Dmol
    _PY3DMOL_AVAILABLE = True
except ImportError:
    py3Dmol = None
    _PY3DMOL_AVAILABLE = False


class MolecularVisualizer:
    """
    Визуализация молекулярных структур.
    
    Методы:
    - show_structure: 3D интерактивная визуализация
    - show_trajectory: Анимация MD-траектории
    - save_image: Сохранение скриншотов
    """
    
    @staticmethod
    def show_structure(
        pdb_file: str,
        style: str = "cartoon",
        color_scheme: str = "spectrum",
        width: int = 800,
        height: int = 600,
        surface_opacity: float = 0.7
    ):
        """
        Интерактивная 3D визуализация структуры.
        
        Args:
            pdb_file: Путь к PDB
            style: Стиль отображения (cartoon, stick, sphere, line)
            color_scheme: Цветовая схема (spectrum, chain, residue, etc.)
            width, height: Размеры окна
            surface_opacity: Прозрачность поверхности (0-1)
        
        Returns:
            py3Dmol view объект (для Jupyter)
        """
        if not _PY3DMOL_AVAILABLE:
            raise ImportError(
                "Для 3D-визуализации установите: pip install py3Dmol"
            )

        with open(pdb_file) as f:
            pdb_string = f.read()
        
        view = py3Dmol.view(width=width, height=height)
        view.addModel(pdb_string, "pdb")
        
        # Стиль для белка
        view.setStyle(
            {"chain": "A"},
            {style: {"color": color_scheme}}
        )
        
        # Стиль для лиганда
        view.setStyle(
            {"chain": "L"},
            {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}}
        )
        
        # Поверхность
        view.addSurface(
            py3Dmol.VDW,
            {"opacity": surface_opacity, "color": "white"},
            {"chain": "A"}
        )
        
        view.zoomTo()
        return view
    
    @staticmethod
    def show_trajectory(
        trajectory_file: str,
        topology_file: str,
        frame_skip: int = 10,
        width: int = 800,
        height: int = 600
    ):
        """
        Анимация траектории (для Jupyter).
        
        Args:
            trajectory_file: Путь к траектории
            topology_file: Путь к топологии
            frame_skip: Показывать каждый N-й кадр
        """
        if not _PY3DMOL_AVAILABLE:
            raise ImportError(
                "Для 3D-визуализации установите: pip install py3Dmol"
            )

        traj = md.load(trajectory_file, top=topology_file)
        
        # Подвыборка кадров
        traj = traj[::frame_skip]
        
        view = py3Dmol.view(width=width, height=height)
        
        for frame in traj:
            pdb_buffer = io.StringIO()
            frame.save_pdb(pdb_buffer)
            view.addModel(pdb_buffer.getvalue(), "pdb")
        
        view.setStyle({"cartoon": {"color": "spectrum"}})
        view.animate({"loop": "forward", "interval": 100})
        view.zoomTo()
        
        return view
    
    @staticmethod
    def create_comparison_plot(
        original_pdb: str,
        simulated_pdb: str,
        output_file: str = "comparison.png"
    ):
        """
        Сравнение начальной и конечной структуры.
        """
        fig = plt.figure(figsize=(15, 7))
        
        # Оригинальная структура
        ax1 = fig.add_subplot(121, projection='3d')
        MolecularVisualizer._plot_structure_3d(original_pdb, ax1, "Начальная структура")
        
        # Симулированная структура
        ax2 = fig.add_subplot(122, projection='3d')
        MolecularVisualizer._plot_structure_3d(simulated_pdb, ax2, "После симуляции")
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"[Viz] Сохранено сравнение: {output_file}")
    
    @staticmethod
    def _plot_structure_3d(pdb_file: str, ax, title: str):
        """Вспомогательная функция для 3D matplotlib."""
        coords = []
        
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM") and " CA " in line:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
        
        coords = np.array(coords)
        
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 
                'o-', markersize=2, linewidth=0.5, alpha=0.6)
        ax.set_title(title)
        ax.set_xlabel("X (Å)")
        ax.set_ylabel("Y (Å)")
        ax.set_zlabel("Z (Å)")
    
    @staticmethod
    def plot_contact_map(
        pdb_file: str,
        cutoff: float = 8.0,
        output_file: str = "contact_map.png"
    ):
        """
        Карта контактов (расстояния между C-alpha).
        """
        traj = md.load(pdb_file)
        
        # Выбор C-alpha
        ca_indices = traj.topology.select("name CA")
        
        # Матрица расстояний C-alpha (нм) - считается напрямую,
        # т.к. md.compute_contacts пропускает нестандартные остатки
        coords = traj.xyz[0, ca_indices]
        dist_matrix = cdist(coords, coords)
        n = len(ca_indices)
        
        # Бинаризация
        contact_matrix = (dist_matrix < cutoff / 10).astype(int)  # Å → нм
        
        # Визуализация
        plt.figure(figsize=(10, 8))
        plt.imshow(contact_matrix, cmap='Blues', interpolation='nearest')
        plt.colorbar(label='Контакт')
        plt.xlabel('Остаток')
        plt.ylabel('Остаток')
        plt.title(f'Карта контактов (cutoff={cutoff} Å)')
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        print(f"[Viz] Карта контактов: {output_file}")


# ============ ТЕСТ ============
if __name__ == "__main__":
    viz = MolecularVisualizer()
    
    # Для Jupyter:
    # view = viz.show_structure("test_data/complex.pdb")
    # view.show()
    
    # Для standalone:
    viz.create_comparison_plot(
        "test_data/initial.pdb",
        "test_data/final.pdb"
    )
