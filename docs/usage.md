# Usage

All pipeline operations go through the `digraph` CLI, which wraps Snakemake and validates your configuration before dispatching any jobs.

---

## Global options

Every command accepts the following options:

| Option | Short | Description |
| --- | --- | --- |
| `--config PATH` | `-c` | Path to `config.yaml` **(required)** |
| `--cores N` | `-t` | Override `resources.snakemake_cores` from the config |
| `--dry-run` | `-n` | Print the execution plan without running any jobs |
| `--force` | `-f` | Force re-execution of all rules regardless of existing outputs |
| `--verbose` | `-v` | Enable DEBUG-level logging |

---

## `digraph validate`

Parse and validate the configuration file. Reports missing FASTQ files and prints a summary of samples, timepoints, and resources. Does not launch any pipeline jobs.

```bash
digraph validate --config config/config.yaml
```

Example output:

```
Config is valid.
  Samples   : 1_Wt, 2_exo1, 3_sgs1
  Timepoints: T0, TSG, TLG, TLR
  Replicates: E1, E2, E3
  Cores     : 8
```

Run this before every new experiment to catch configuration or data-preparation errors early.

---

## `digraph run`

Run the complete pipeline from start to finish.

```bash
# Full run on 8 cores
digraph run --config config/config.yaml --cores 8

# Preview the DAG without executing anything
digraph run --config config/config.yaml --dry-run

# Force all rules to re-run (ignore existing outputs)
digraph run --config config/config.yaml --cores 8 --force

# Stop at a specific Snakemake rule
digraph run --config config/config.yaml --until r_plot_gal_vs_raf
```

Snakemake automatically resumes from where it stopped if the run is interrupted — no jobs that already completed will re-run.

---

## `digraph stage`

Run the pipeline up to and including a single named stage. Useful for incremental analysis or debugging.

```bash
digraph stage coverage    --config config/config.yaml --cores 8
digraph stage categories  --config config/config.yaml --cores 8
digraph stage mutagenic   --config config/config.yaml --cores 8
digraph stage discordant  --config config/config.yaml --cores 8
digraph stage report      --config config/config.yaml --cores 8
```

Each stage command is equivalent to `digraph run --until <terminal_rule>`:

| Stage name | Terminal Snakemake rule |
| --- | --- |
| `coverage` | `r_process_cov_18nt` |
| `categories` | `r_plot_gal_vs_raf` |
| `mutagenic` | `r_plot_repair_comparison` |
| `discordant` | `r_discordant_network` |
| `report` | `generate_report` |

---

## Cluster execution

Di-GRAPH runs locally by default. To dispatch jobs to an HPC cluster, install a [Snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/) and pass extra arguments via `digraph run`:

```bash
# SLURM example
digraph run --config config/config.yaml \
    -- --executor slurm --jobs 32 --default-resources mem_mb=8000 runtime=120
```

Arguments after `--` are forwarded verbatim to Snakemake.

---

## Typical workflow

```bash
# 1. Validate config and check inputs
digraph validate --config config/config.yaml

# 2. Preview the full DAG
digraph run --config config/config.yaml --dry-run

# 3. Run coverage stage first (fastest, good sanity check)
digraph stage coverage --config config/config.yaml --cores 8

# 4. Run the rest
digraph run --config config/config.yaml --cores 8
```
