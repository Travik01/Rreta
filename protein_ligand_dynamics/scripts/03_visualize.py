#!/usr/bin/env python3
"""
Точка входа визуализации (совместимость с исходной нумерацией файлов).

Использование:
    python scripts/03_visualize.py --pdb results/complex.pdb --trajectory results/trajectory.dcd
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from visualize import main  # noqa: E402

if __name__ == "__main__":
    main()
