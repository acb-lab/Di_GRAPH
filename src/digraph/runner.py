"""
Programmatic Snakemake invocation for the Di-GRAPH pipeline.

Builds the ``snakemake`` command from a validated :class:`~digraph.config.DiGraphConfig`
and runs it via :func:`subprocess.run`.  Using subprocess rather than the
Snakemake Python API keeps the CLI package decoupled from any specific Snakemake
version — only the installed ``snakemake`` binary on ``PATH`` matters.

The module also exposes :data:`STAGE_TERMINAL_RULES`, which maps human-readable
stage names to the last Snakemake rule in each stage.  This is used by
``digraph stage`` to pass ``--until <rule>`` to Snakemake.

Usage
-----
::

    from pathlib import Path
    from digraph.config import load_config
    from digraph.runner import run_snakemake

    config = load_config(Path("config/config.yaml"))
    exit_code = run_snakemake(
        config,
        Path("config/config.yaml"),
        cores=8,
        dry_run=False,
    )
"""

from __future__ import annotations

import json
import logging
import subprocess
import tempfile
from pathlib import Path

from digraph.config import DiGraphConfig
from digraph.utils.io import snakemake_config_dict

__all__ = ["STAGE_TERMINAL_RULES", "run_snakemake"]

logger = logging.getLogger("digraph")

#: Maps human-readable stage names to the terminal Snakemake rule of that stage.
#: Passed to ``snakemake --until <rule>`` by :func:`run_snakemake` when a
#: specific stage is requested via ``digraph stage <name>``.
STAGE_TERMINAL_RULES: dict[str, str] = {
    "coverage":   "r_process_cov_18nt",       # Stage 1 — coverage tracks + MAT quant
    "categories": "r_plot_gal_vs_raf",         # Stage 2 — genomic category fingerprints
    "mutagenic":  "r_plot_repair_comparison",  # Stage 3 — repair pathway choice
    "discordant": "r_discordant_network",      # Stage 4 — genome-wide rearrangements
    "report":     "generate_report",           # Stage 5 — HTML report
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
