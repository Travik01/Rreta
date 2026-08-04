#!/usr/bin/env python3
"""
Точка входа пайплайна (совместимость с исходной нумерацией файлов).

Использование:
    python scripts/01_full_pipeline.py --config config/pipeline_config.yaml
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from full_pipeline import main  # noqa: E402

if __name__ == "__main__":
    main()
