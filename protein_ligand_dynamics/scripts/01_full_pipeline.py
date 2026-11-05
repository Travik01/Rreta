#!/usr/bin/env python3
# ruff: noqa: N999
"""Точка входа полного пайплайна."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from full_pipeline import main

if __name__ == "__main__":
    main()
