"""
Path helpers and input-file validation for Di-GRAPH.

This module provides three public utilities consumed by the CLI before any
Snakemake job is dispatched:

:func:`check_fastq_inputs`
    Enumerate every expected paired FASTQ.gz file and return a list of
    those that are missing.  Called by ``digraph run`` and ``digraph validate``
    to surface data-preparation errors before the pipeline starts.

:func:`build_sample_map`
    Build a ``{name: path}`` dict for quick lookups of sample directories.

:func:`snakemake_config_dict`
    Serialise the validated :class:`~digraph.config.DiGraphConfig` to a plain
    ``dict`` (converting :class:`~pathlib.Path` objects to strings) suitable
    for passing to Snakemake as a ``--config`` value or a temporary config file.
"""

from __future__ import annotations

import logging
from pathlib import Path

from digraph.config import DiGraphConfig, ExperimentConfig, SampleConfig

__all__ = ["check_fastq_inputs", "build_sample_map", "snakemake_config_dict"]

logger = logging.getLogger("digraph")


def check_fastq_inputs(config: DiGraphConfig) -> list[str]:
    """
    Verify that all expected paired FASTQ.gz files exist before launching
    Snakemake.  Returns a list of missing-file messages; empty on success.

    Expected file name convention::

        <working_dir>/<strain>/<timepoint>_<experiment>_R1.fastq.gz
        <working_dir>/<strain>/<timepoint>_<experiment>_R2.fastq.gz

    Args:
        config: Validated :class:`~digraph.config.DiGraphConfig` instance.

    Returns:
        List of human-readable error strings; callers should raise if non-empty.
    """
    missing: list[str] = []
    for sample in config.samples:
        for timepoint in config.experiments.timepoints:
            for exp in config.experiments.names:
                for read in ("R1", "R2"):
                    fq = sample.path / f"{timepoint}_{exp}_{read}.fastq.gz"
                    if not fq.exists():
                        missing.append(str(fq))
    return missing


def build_sample_map(config: DiGraphConfig) -> dict[str, Path]:
    """
    Build a mapping of sample name → sample directory path.

    Args:
        config: Validated :class:`~digraph.config.DiGraphConfig`.

    Returns:
        Dict such as ``{"1_Wt": Path("/data/.../1_Wt"), ...}``.
    """
    return {s.name: s.path for s in config.samples}


def snakemake_config_dict(config: DiGraphConfig) -> dict:
    """
    Flatten the validated config into a plain dict suitable for passing to
    Snakemake via ``--config`` or a temporary config file.

    All ``Path`` objects are converted to strings so Snakemake can serialise
    them as JSON / YAML.

    Args:
        config: Validated :class:`~digraph.config.DiGraphConfig`.

    Returns:
        Nested dict mirroring the config structure with paths as strings.
    """
    return _paths_to_str(config.model_dump())


def _paths_to_str(obj: object) -> object:
    """Recursively convert Path objects to strings inside nested dicts/lists."""
    if isinstance(obj, dict):
        return {k: _paths_to_str(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_paths_to_str(item) for item in obj]
    if isinstance(obj, Path):
        return str(obj)
    return obj
