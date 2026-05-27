"""
Di-GRAPH — DSB-induced Genome-wide Repair Analysis and Profiling of Homologous Recombination.

This package provides the ``digraph`` command-line interface and a programmatic API
for running the five-stage Di-GRAPH bioinformatics pipeline.

Pipeline stages
---------------
1. **Coverage** — Read trimming, alignment (bowtie), coverage tracks (deeptools),
   and MAT-locus polymorphism quantification.
2. **Categories** — Per-category coverage fingerprints across 13 genomic feature
   classes (ORF, LTR, TEG, Ty, tRNA, rRNA, ncRNA, snRNA, snoRNA, ARS, Cen, Tel, Int).
3. **Mutagenic** — BWA alignment, HO-site read extraction, and repair pathway
   classification (gene conversion vs NHEJ vs other).
4. **Discordant** — Inter-chromosomal discordant read mapping, BLAST cross-validation,
   and genome-wide rearrangement profiling.
5. **Report** — RMarkdown flexdashboard HTML report summarising all stages.

Package layout
--------------
``digraph.config``
    Pydantic v2 models for ``config.yaml`` validation and loading.
    Entry point: :func:`~digraph.config.load_config`.

``digraph.cli``
    Typer-based CLI — ``digraph run``, ``digraph validate``, ``digraph stage``.

``digraph.runner``
    Subprocess wrapper that builds and executes the ``snakemake`` command.
    Entry point: :func:`~digraph.runner.run_snakemake`.

``digraph.utils.logging``
    Rich-based logging initialisation.
    Entry point: :func:`~digraph.utils.logging.setup_logging`.

``digraph.utils.io``
    FASTQ input validation and config serialisation helpers.

Quick start (CLI)
-----------------
After ``pip install -e .`` inside the activated ``digraph`` conda environment::

    # Validate config without running anything
    digraph validate --config config/config.yaml

    # Preview the full execution DAG
    digraph run --config config/config.yaml --dry-run

    # Run the complete pipeline on 8 cores
    digraph run --config config/config.yaml --cores 8

    # Run only stage 1 (coverage)
    digraph stage coverage --config config/config.yaml --cores 8

Programmatic usage
------------------
::

    from pathlib import Path
    from digraph.config import load_config
    from digraph.runner import run_snakemake

    config = load_config(Path("config/config.yaml"))
    exit_code = run_snakemake(config, Path("config/config.yaml"), dry_run=True)
"""

__version__ = "0.5.0"

__all__ = ["__version__"]
