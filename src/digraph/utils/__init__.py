"""
digraph.utils — shared utility helpers.

This sub-package collects small, focused modules that are used by both the CLI
and the Snakemake runner.  None of the modules here contain pipeline logic.

Modules
-------
``digraph.utils.logging``
    Initialise a ``rich``-based logger for the ``digraph`` namespace.
    Used by all CLI commands to provide coloured, structured console output
    and optional file logging.

``digraph.utils.io``
    Pre-flight FASTQ validation, sample-map construction, and config
    serialisation to a plain dict for Snakemake.
"""
