"""Настройка логирования."""

import logging
import sys
from pathlib import Path


def setup_logger(
    name: str,
    log_file: str | Path | None = None,
    level: int = logging.INFO,
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    has_console = any(
        isinstance(handler, logging.StreamHandler)
        and not isinstance(handler, logging.FileHandler)
        for handler in logger.handlers
    )
    if not has_console:
        console = logging.StreamHandler(sys.stdout)
        console.setFormatter(formatter)
        logger.addHandler(console)

    if log_file:
        path = Path(log_file).resolve()
        has_file = any(
            isinstance(handler, logging.FileHandler)
            and Path(handler.baseFilename).resolve() == path
            for handler in logger.handlers
        )
        if not has_file:
            path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(str(path), encoding="utf-8")
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

    return logger
