"""Logging configuration for MPLID pipeline.

Author: Folorunsho Bright Omage
License: MIT
"""

import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional


def setup_pipeline_logging(
    name: str,
    log_dir: Optional[Path] = None,
    level: int = logging.INFO,
    console: bool = True,
) -> logging.Logger:
    """Configure logging for pipeline scripts.

    Parameters
    ----------
    name : str
        Logger name (typically __name__).
    log_dir : Optional[Path], optional
        Directory for log files. If None, file logging is disabled.
    level : int, optional
        Logging level (default: logging.INFO).
    console : bool, optional
        Whether to log to console (default: True).

    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Clear existing handlers
    logger.handlers.clear()

    # Formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Console handler
    if console:
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # File handler
    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = log_dir / f"{name}_{timestamp}.log"

        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


class PipelineProgress:
    """Simple progress tracking for pipeline stages."""

    def __init__(self, logger: logging.Logger, total: int, desc: str = "Progress"):
        """Initialize progress tracker.

        Parameters
        ----------
        logger : logging.Logger
            Logger instance for output.
        total : int
            Total number of items to process.
        desc : str, optional
            Description of the progress (default: "Progress").
        """
        self.logger = logger
        self.total = total
        self.desc = desc
        self.current = 0
        self.last_logged = 0

    def update(self, n: int = 1) -> None:
        """Update progress by n items.

        Parameters
        ----------
        n : int, optional
            Number of items completed (default: 1).
        """
        self.current += n
        percent = (self.current / self.total) * 100

        # Log at 10% intervals
        if int(percent / 10) > int(self.last_logged / 10):
            self.logger.info(f"{self.desc}: {self.current}/{self.total} ({percent:.1f}%)")
            self.last_logged = percent

    def finish(self) -> None:
        """Mark progress as complete."""
        self.logger.info(f"{self.desc}: Complete ({self.total}/{self.total})")
