# Installation

Di-GRAPH requires **conda** to manage the bioinformatics runtime (R 4.4, Python 3.9, bowtie, samtools, deeptools, BLAST, etc.) and **pip** to install the `digraph` CLI into that environment.

---

## Requirements

| Requirement | Version | Purpose |
| --- | --- | --- |
| conda / mamba | any | environment manager |
| Python | ≥ 3.11 | `digraph` CLI |
| Git | any | cloning the repository |

All bioinformatics tools (bowtie, bowtie2, bwa, samtools, deeptools, bcftools, BLAST, R, pysam, …) are installed automatically by the conda environment file and do not need to be installed manually.

---

## Step 1 — Clone the repository

```bash
git clone https://github.com/acb-lab/Di_GRAPH.git
cd Di_GRAPH
```

---

## Step 2 — Create the conda environment

```bash
conda env create -f environment/digraph.yml -n digraph
conda activate digraph
```

This installs all tools and R/Python packages listed in `environment/digraph.yml`, including R 4.4 with tidyverse, ggplot2, ggraph, flexdashboard, and all other packages required by the analysis scripts.

!!! tip
    Using [mamba](https://mamba.readthedocs.io/) instead of conda significantly speeds up environment creation:
    ```bash
    mamba env create -f environment/digraph.yml -n digraph
    ```

---

## Step 3 — Install the `digraph` CLI

```bash
pip install -e .
```

This installs the `digraph` command into the active conda environment. The `-e` flag makes it an editable install, so local changes to `src/digraph/` take effect immediately without reinstalling.

Verify the installation:

```bash
digraph --help
```

---

## Developer installation

For contributors who also want linting and type-checking tools:

```bash
pip install -e ".[dev]"
```

For building the documentation locally:

```bash
pip install -e ".[docs]"
mkdocs serve
```

Then open <http://127.0.0.1:8000> in your browser.

---

## Updating

To update Di-GRAPH after pulling a new version:

```bash
git pull
conda env update -f environment/digraph.yml -n digraph --prune
pip install -e .
```
