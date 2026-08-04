#!/usr/bin/env python3
"""
Точка входа анализа (совместимость с исходной нумерацией файлов).

Использование:
    python scripts/02_analyze_results.py --trajectory results/trajectory.dcd \
                                         --topology results/complex.pdb
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from analyze_results import main  # noqa: E402

if __name__ == "__main__":
    main()
