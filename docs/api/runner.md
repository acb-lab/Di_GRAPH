# Runner

Subprocess-based Snakemake invocation. Builds the `snakemake` command from a validated config and runs it via `subprocess.run`.

---

::: digraph.runner
    options:
      members:
        - STAGE_TERMINAL_RULES
        - run_snakemake
      show_source: false
      show_root_heading: true
