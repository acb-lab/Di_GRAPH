"""Structured logging setup using the ``rich`` library."""

import logging
from pathlib import Path

from rich.logging import RichHandler


def setup_logging(log_file: Path | None = None, verbose: bool = False) -> logging.Logger:
    """
    Configure the root logger with a Rich console handler and an optional
    file handler.

    Args:
        log_file: If provided, log messages are also written to this file.
        verbose:  When True, sets level to DEBUG; otherwise INFO.

    Returns:
        The ``digraph`` logger instance.
    """
    level = logging.DEBUG if verbose else logging.INFO

    handlers: list[logging.Handler] = [
        RichHandler(rich_tracebacks=True, markup=True, show_path=False)
    ]

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        )
        handlers.append(file_handler)

    logging.basicConfig(level=level, handlers=handlers, force=True)
    return logging.getLogger("digraph")
