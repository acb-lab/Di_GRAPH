"""
Di-GRAPH command-line interface.

Provides three commands:

* ``digraph run``      — run the full pipeline
* ``digraph validate`` — validate the config file without running anything
* ``digraph stage``    — run a single named pipeline stage

Usage (after ``pip install -e .`` or ``uv pip install -e .``):

    digraph run --config config/config.yaml --cores 8
    digraph validate --config config/config.yaml
    digraph stage coverage --config config/config.yaml --dry-run
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console

from digraph.config import load_config
from digraph.runner import STAGE_TERMINAL_RULES, run_snakemake
from digraph.utils.io import check_fastq_inputs
from digraph.utils.logging import setup_logging

app = typer.Typer(
    name="digraph",
    help="Di-GRAPH: DSB-induced Genome-wide Repair Analysis pipeline.",
    add_completion=False,
)

console = Console()

# ---------------------------------------------------------------------------
# Shared option types
# ---------------------------------------------------------------------------

_ConfigOption = typer.Option(..., "--config", "-c", help="Path to config.yaml")
_CoresOption = typer.Option(None, "--cores", "-t", help="Override snakemake_cores from config")
_DryRunOption = typer.Option(False, "--dry-run", "-n", help="Print rules without executing")
_VerboseOption = typer.Option(False, "--verbose", "-v", help="Enable debug logging")


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


@app.command()
def run(
    config_path: Path = _ConfigOption,
    cores: Optional[int] = _CoresOption,
    dry_run: bool = _DryRunOption,
    force: bool = typer.Option(False, "--force", "-f", help="Force re-execution of all rules"),
    until: Optional[str] = typer.Option(None, "--until", help="Run up to this Snakemake rule"),
    verbose: bool = _VerboseOption,
) -> None:
    """Run the full Di-GRAPH pipeline."""
    log = setup_logging(verbose=verbose)

    config = _load_or_exit(config_path, log)
    _check_inputs_or_exit(config, log)

    code = run_snakemake(
        config,
        config_path,
        dry_run=dry_run,
        force_rerun=force,
        until_rule=until,
        cores=cores,
    )
    sys.exit(code)


@app.command()
def validate(
    config_path: Path = _ConfigOption,
    verbose: bool = _VerboseOption,
) -> None:
    """Validate the config.yaml file without running the pipeline."""
    log = setup_logging(verbose=verbose)
    config = _load_or_exit(config_path, log)

    missing = check_fastq_inputs(config)
    if missing:
        log.warning("The following expected FASTQ files are missing:")
        for f in missing:
            log.warning("  %s", f)
    else:
        log.info("All expected FASTQ inputs found.")

    console.print("[bold green]Config is valid.[/bold green]")
    console.print(f"  Samples   : {', '.join(config.sample_names)}")
    console.print(f"  Timepoints: {', '.join(config.experiments.timepoints)}")
    console.print(f"  Replicates: {', '.join(config.experiments.names)}")
    console.print(f"  Cores     : {config.resources.snakemake_cores}")


@app.command()
def stage(
    stage_name: str = typer.Argument(
        help=f"Stage to run: {', '.join(STAGE_TERMINAL_RULES)}"
    ),
    config_path: Path = _ConfigOption,
    cores: Optional[int] = _CoresOption,
    dry_run: bool = _DryRunOption,
    force: bool = typer.Option(False, "--force", "-f"),
    verbose: bool = _VerboseOption,
) -> None:
    """Run a single named pipeline stage."""
    log = setup_logging(verbose=verbose)

    if stage_name not in STAGE_TERMINAL_RULES:
        valid = ", ".join(STAGE_TERMINAL_RULES)
        log.error("Unknown stage '%s'. Valid stages: %s", stage_name, valid)
        sys.exit(1)

    config = _load_or_exit(config_path, log)
    _check_inputs_or_exit(config, log)

    until_rule = STAGE_TERMINAL_RULES[stage_name]
    log.info("Running stage '%s' (until rule: %s)", stage_name, until_rule)

    code = run_snakemake(
        config,
        config_path,
        dry_run=dry_run,
        force_rerun=force,
        until_rule=until_rule,
        cores=cores,
    )
    sys.exit(code)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_or_exit(config_path: Path, log) -> object:  # type: ignore[return]
    """Load config, printing a friendly error and exiting on failure."""
    try:
        from digraph.config import DiGraphConfig  # local import for clarity

        cfg = load_config(config_path)
        log.info("Config loaded: %s", config_path)
        return cfg
    except Exception as exc:
        console.print(f"[bold red]Config error:[/bold red] {exc}")
        sys.exit(1)


def _check_inputs_or_exit(config, log) -> None:
    """Warn about missing FASTQ inputs; does not abort (some may be intentional)."""
    missing = check_fastq_inputs(config)
    if missing:
        log.warning("%d expected FASTQ files not found (first 5 shown):", len(missing))
        for f in missing[:5]:
            log.warning("  %s", f)
