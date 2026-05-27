"""
Programmatic Snakemake invocation.

Builds the ``snakemake`` command from a validated :class:`DiGraphConfig`
and runs it via :func:`subprocess.run`.  Using subprocess rather than the
Snakemake Python API avoids tight version coupling between the CLI package
and the Snakemake library.
"""

from __future__ import annotations

import json
import logging
import subprocess
import tempfile
from pathlib import Path

from digraph.config import DiGraphConfig
from digraph.utils.io import snakemake_config_dict

logger = logging.getLogger("digraph")

# Map stage names to the terminal Snakemake rule for ``--until``
STAGE_TERMINAL_RULES: dict[str, str] = {
    "coverage": "r_process_cov_18nt",
    "categories": "r_plot_gal_vs_raf",
    "mutagenic": "r_plot_repair_comparison",
    "discordant": "r_discordant_network",
    "report": "generate_report",
}

# Path to the top-level Snakefile relative to this file's package root
_PACKAGE_ROOT = Path(__file__).parent.parent.parent
_SNAKEFILE = _PACKAGE_ROOT / "workflow" / "Snakefile"


def run_snakemake(
    config: DiGraphConfig,
    config_path: Path,
    *,
    dry_run: bool = False,
    force_rerun: bool = False,
    until_rule: str | None = None,
    cores: int | None = None,
    extra_args: list[str] | None = None,
) -> int:
    """
    Invoke Snakemake as a subprocess with the given configuration.

    Args:
        config:       Validated pipeline configuration.
        config_path:  Path to the config YAML file (passed as ``--configfile``).
        dry_run:      Pass ``--dryrun`` to Snakemake.
        force_rerun:  Pass ``--forceall`` to Snakemake.
        until_rule:   Stop after this rule (``--until``).
        cores:        Override ``config.resources.snakemake_cores``.
        extra_args:   Additional raw arguments appended to the command.

    Returns:
        The process return code (0 = success).
    """
    n_cores = cores if cores is not None else config.resources.snakemake_cores

    cmd: list[str] = [
        "snakemake",
        "--snakefile", str(_SNAKEFILE),
        "--cores", str(n_cores),
        "--configfile", str(config_path),
        "--rerun-incomplete",
        "--printshellcmds",
    ]

    if config.resources.use_conda:
        cmd.append("--use-conda")

    if dry_run:
        cmd.append("--dryrun")

    if force_rerun:
        cmd.append("--forceall")

    if until_rule:
        cmd += ["--until", until_rule]

    if extra_args:
        cmd += extra_args

    logger.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd)
    return result.returncode
