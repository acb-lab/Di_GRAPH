# Di-GRAPH

``` text
██████╗               ██████╗ ██████╗  █████╗ ██████╗ ██╗  ██╗        
██╔══██╗  ██╗        ██╔════╝ ██╔══██╗██╔══██╗██╔══██╗██║  ██║
██║  ██║  ╚═╝  ███╗  ██║  ███╗██████╔╝███████║██████╔╝███████║
██║  ██║  ██╗  ╚══╝  ██║   ██║██╔╚██╗ ██╔══██║██╔═══╝ ██╔══██║
██████╔╝  ██║        ╚██████╔╝██║ ╚██╗██║  ██║██║     ██║  ██║
╚═════╝   ╚═╝         ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝

        ===========================================
          DSB-induced Genome-wide Repair Analysis
         and Profiling of Homologous recombination
        =========================================== 
```

<br>

**Di-GRAPH** (<ins>D</ins>SB-<ins>i</ins>nduced <ins>G</ins>enome-wide <ins>R</ins>epair <ins>A</ins>nalysis and <ins>P</ins>rofiling of <ins>H</ins>omologous recombination)

A pipeline to detect, classify and interpret recombination events at a defined break site and across the entire genome upon HO-induction of a single site-specific DSB in *S. cerevisiae*. 
________________________________________________________________________________________________________________________________________________

## Table of contents<a name="idindex"></a>

1.  [Introduction](#idintro)
2.  [Installation](#idinstall)
3.  [Instructions](#idinstr)
4.  [Example usage](#idexample)
5.  [Expected output](#idoutput)
6.  [Optimisation opportunities](#idoptimise)
7.  [Parallelisation](#idparallel)
8.  [Visual summary](#idsummary)

<br>

## 1. Introduction<a name="idintro"></a>

**Di-GRAPH (**<ins>D</ins>SB-<ins>i</ins>nduced <ins>G</ins>enome-wide <ins>R</ins>epair <ins>A</ins>nalysis and <ins>P</ins>rofiling of <ins>H</ins>omologous recombination) integrates DSB mutational signature analysis, repair pathway choice, coverage profiling and discordant read mapping to quantitatively define the frequency, directionality, extent and mutagenic profile of gene conversion events during the repair of a single site-specific DSB.

At a genomic scale, Di-GRAPH identifies and maps DNA damage-dependent genome-wide gross chromosomal rearrangements in defined genomic regions to evaluate how the genome is reshaped in response to a DSB.

This pipeline is designed to analyze how stage-specific HR repair factors differentially control DSB repair fidelity and genome-wide stability in the *S. cerevisiae* PMV genetic background, which allows the induction of a DSB in the *MATa* locus on chromosome III in a galactose-dependent manner; and contains an engineered *MATa'* locus on chromosome V used as donor for recombination. 

> Note: detailed information about the PMV genetic background is available in *Ramos et al.,2022 - Cell Reports*, <https://doi.org/10.1016/j.celrep.2021.110201>


By comparing data prior to DSB induction, during its repair (non-selected survivors) and from survivor populations (selected survivors), Di-GRAPH enables the classification of lethal *vs* non-lethal rearrangements arising during the repair of the DSB. Additionally, Di-GRAPH incorporates the assessment of undamaged contitions (undamaged cells) to distinguish DNA damage-dependent from cell-cycle-dependent genomic alterations.

After providing paired-end genomic sequencing data from these 4 different timepoints, Di-GRAPH will:

- i)  Perform *MATa/MATa'* loci coverage profiling to identify gene conversion patterns and polymporphisms incorporation.

- ii) Characterize the frequency, directionality and extent of individual gene conversion events between *MATa/MATa'* loci by applying inter-chromosomal discordant read mapping.

- iii) Assess HO associated mutagenic pattern and repair pathway choice upon DSB induction in the *MATa* locus.

- iv) Evaluate global genome stability by characterizing how different genomic categories (e.g. ORFs, intergenic regions, LTR, TEG, Ty, tRNA, rRNA, ncRNA, snRNA, snoRNA, ARS, centromere and subtelomeric regions) behave in cells lacking defined HR factors, both after DSB induction and in the absence of DNA damage.

- v) Apply inter-chromosomal discordant read mapping combined with BLAST cross-validation to identify and map genome-wide chromosomal rearrangements in both DNA-damaged and undamaged cells, characterizing how the genome is reshaped in the absence of defined HR factors.

<br>

You can find a **visual summary** of the different steps conducted by Di-GRAPH [here](#idsummary).

<br>

[Back to index](#idindex)

<br>

## 2. Installation<a name="idinstall"></a>

Di-GRAPH requires [conda](https://www.anaconda.com/docs/getting-started/miniconda/install) to manage the bioinformatics runtime environment (R, Python, bowtie, samtools, etc.).

### 1. Clone the repository

```bash
git clone https://github.com/acb-lab/Di_GRAPH.git
cd Di_GRAPH
```

### 2. Create the conda environment

```bash
conda env create -f environment/digraph.yml -n digraph
conda activate digraph
```

### 3. Install the `digraph` CLI

```bash
pip install -e .
```

After these three steps, the `digraph` command is available in your activated conda environment.

<br>

[Back to index](#idindex)

<br>

## 3. Instructions<a name="idinstr"></a>

Di-GRAPH is driven by a YAML configuration file (`config/config.yaml`) and invoked through the `digraph` CLI. As noted in the [Introduction](#idintro), Di-GRAPH requires paired-end genomic data (FASTQ files) from 4 different timepoints. Paired-end reads should be a minimum of 150 bp to support the mutagenic pattern analysis. Three independent replicates per timepoint are required.

All reference files (PMV reference genome, genomic categories annotation, BLAST database, and HTML report template) are available in the `files/` folder of this repository and are referenced by default in `config/config.yaml`.

### Configuring your run

Copy and edit `config/config.yaml`, setting at minimum the `paths.working_dir` and the `samples` list:

```yaml
paths:
  working_dir: /path/to/your/experiment_directory

samples:
  - name: Wt
  - name: exo1
  - name: sgs1
```

All other paths in the config file are relative to the repository root and will work without modification if you run `digraph` from the repository directory. Adjust `resources.snakemake_cores` to match your machine.

### Available CLI commands

```
digraph run      --config config/config.yaml [--cores N] [--dry-run] [--until RULE] [--force]
digraph validate --config config/config.yaml
digraph stage    coverage|categories|mutagenic|discordant|report --config config/config.yaml
```

`digraph validate` checks that all paths exist and the config is complete before any computation starts.  
`digraph stage` runs the pipeline only up to the specified stage.  
`digraph run --dry-run` prints the full execution plan without running any jobs.

<br>

### Input data structure

For each strain listed in the config, a subfolder with the same name must exist inside the working directory. Paired-end FASTQ files go inside that subfolder and must follow the naming convention `timepoint_replicate_R1.fastq.gz` / `timepoint_replicate_R2.fastq.gz`:

```text
------- Working directory (created by the user) -------
working_dir/
├── Wt/
│   ├── T0_E1_R1.fastq.gz
│   ├── T0_E1_R2.fastq.gz
│   ├── TSG_E1_R1.fastq.gz
│   ├── TSG_E1_R2.fastq.gz
│   └── ...  (TLG_E1, TLR_E1, T0_E2, TSG_E2, ...)
├── exo1/
│   └── ...
└── Strain_n/
    └── ...

------- Reference files (available in Di_GRAPH/files/) -------
files/
├── RG/                          # PMV reference genome
│   ├── RG_PMV_v9.fasta
│   └── ...
├── Categories/                  # genomic categories annotation
│   ├── PMV_categories.tsv
│   └── ...
├── BLAST/features_extraction/   # BLAST cross-validation database
│   ├── YLL039C_seq.fasta
│   └── ...
└── Report_files/                # HTML report template
    ├── Di-GRAPH_report.Rmd
    └── ...
```

> **Timepoint names:** **T0** (prior to DSB induction), **TSG** (Short Galactose — non-selected survivors), **TLG** (Long Galactose — selected survivors), **TLR** (Long Raffinose — undamaged cells).  
> **Replicate names:** **E1**, **E2**, **E3**.

<br>

You can see an example of Di-GRAPH usage in [this section](#idexample).

<br>

[Back to index](#idindex)

<br>

## 4. Example usage<a name="idexample"></a>

You can test Di-GRAPH with paired-end genomic data from wild-type, *exo1∆*, *sgs1∆*, *srs2∆* and *rad51∆* PMV cells included in the `test_dataset` directory.

**1. Edit the config** to point `working_dir` at the test dataset and list the sample names:

```yaml
# config/config.yaml
paths:
  working_dir: test_dataset/working_directory

samples:
  - name: 1_Wt
    is_reference: true
  - name: 2_exo1
  - name: 3_sgs1
  - name: 4_srs2
  - name: 5_rad51
```

### 2. Validate the configuration

```bash
digraph validate --config config/config.yaml
```

### 3. Preview the execution plan (dry run)

```bash
digraph run --config config/config.yaml --dry-run
```

### 4. Run the full pipeline

```bash
digraph run --config config/config.yaml --cores 8
```

To run only a specific stage (e.g. coverage and categories):

```bash
digraph stage coverage   --config config/config.yaml --cores 8
digraph stage categories --config config/config.yaml --cores 8
```

<br>


[Back to index](#idindex)

<br>

## 5. Expected output<a name="idoutput"></a>

For each strain defined in the working directory, Di-GRAPH will perform `bowtie/bowtie2/bwa` genomic alignments, characterize gene conversion products between *MATa/MATa'* loci, define HO associated mutagenic pattern, analyze coverage data regarding all genomic categories and identify global genomic rearrangements. The output files and plots will be stored in the `Working_directory/Strain_n` subfolder(s).

To facilitate the interpretation of the results, Di-GRAPH will generate a `Di-GRAPH_report.html` file in the `Working_directory` folder. This report is a summary of the results obtained by Di-GRAPH for each strain and includes the following sections:
- **Overview:** contains the alignment statistics and the script log file.
- **MAT analysis:** contains coverage analysis and gene conversion analysis for the *MATa/MATa'* loci. It also contains the polymorphisms incorporation analysis.
- **Mutagenic profiling at HO site:** contains the HO associated mutagenic pattern and repair pathway choice analysis.
- **Genome-wide analysis** contains the coverage analysis for all genomic categories and the inter-chromosomal discordant read mapping analysis.

<br>

The Di-GRAPH report file generated after the analysis of paired-end genomic data from wild-type, *exo1∆*, *sgs1∆*, *srs2∆* and *rad51∆* PMV cells is included in the `test_dataset` directory. 

<br>

[Back to index](#idindex)

<br>

## 6. Optimisation opportunities<a name="idoptimise"></a>

The following areas have been identified as candidates for performance improvement in future development iterations.

### Read trimming chain (Stage 1)

The 75 nt coverage analysis requires extracting eight 75 bp fragments from each 150 bp paired-end read. Currently this is done with eight sequential `cutadapt` calls whose outputs are concatenated. A single `fastp` call with a sliding-window extraction parameter would reduce I/O and startup overhead by roughly 8×, and would also allow per-fragment quality trimming to be applied in one pass.

### R script startup overhead (Stages 2–4)

Each R rule spawns an independent R session that reloads the full package stack (tidyverse, ggplot2, ggraph, etc.) from scratch. On a cold filesystem this can dominate the runtime of short-running rules. Possible mitigations:

- Bundle logically related scripts into fewer, larger R jobs.
- Use [`callr`](https://callr.r-lib.org/) to reuse a persistent R subprocess across related rules within a stage.

### BLAST cross-validation (Stage 4)

The BLAST cross-validation iterates over all per-feature FASTA files inside a `while read` shell loop, running one `blastn` call at a time. Because each feature is independent, splitting this loop into individual Snakemake rules (one per FASTA file) would expose full parallelism at no algorithmic cost and is the highest-impact single change available in Stage 4.

### Concordant/discordant SAM processing (Stage 4)

The extraction of inter-discordant read pairs currently relies on a multi-stage `awk` pipeline operating on SAM text. Replacing this with a `pysam`-based Python script would give binary BAM I/O, eliminate intermediate text conversion, and enable per-read filtering logic to be unit-tested.

### Polymorphism coverage calculation (Stage 1)

The per-position baseline subtraction is already implemented in Python (replacing the original `awk` loop), but it reads coverage bedGraph files line-by-line. Switching to `pandas` vectorised operations would reduce memory allocations and improve speed on large coverage files.

<br>

[Back to index](#idindex)

<br>

## 7. Parallelisation<a name="idparallel"></a>

### What runs in parallel today

Di-GRAPH uses [Snakemake](https://snakemake.readthedocs.io/) to manage execution. Snakemake builds a directed acyclic graph (DAG) of all jobs and automatically dispatches every job whose dependencies are satisfied, up to the `--cores` limit. No manual coordination is needed.

With the default dataset of 5 strains × 4 timepoints × 3 replicates, the following jobs are fully independent and run concurrently:

| Stage | Independent unit | Concurrent jobs (5 strains) |
| --- | --- | --- |
| 1 — Coverage | strain × timepoint × replicate | up to 60 alignments |
| 2 — Categories | strain × category (13) | up to 65 fingerprint jobs |
| 3 — Mutagenic | strain × timepoint (T0/TLG/TLR) × replicate | up to 45 BWA jobs |
| 4 — Discordant | strain × timepoint × replicate | up to 60 bowtie2 jobs |

To take full advantage of this, set `resources.snakemake_cores` in `config/config.yaml` to match the number of cores available on your machine, or pass `--cores N` at runtime:

```bash
digraph run --config config/config.yaml --cores 16
```

Individual tools (bowtie, bowtie2, bwa, fastp, bamCoverage) also use multi-threading internally. Their thread count is controlled by `resources.threads` in the config, independently of the Snakemake-level parallelism.

### Cluster and cloud execution

Snakemake supports submitting each rule as an independent job to an HPC scheduler or a cloud provider with no changes to the workflow rules. This is the most impactful scaling option for large sample sets:

- **SLURM / SGE / PBS** — install the corresponding [Snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/) and add a `--executor` flag. Each rule becomes a cluster job with its own resource request.
- **Cloud (AWS, GCP, Azure)** — Snakemake cloud executor plugins allow transparent job dispatch and output storage in object buckets.

```bash
# Example: SLURM cluster with 32 simultaneous jobs
digraph run --config config/config.yaml \
    --executor slurm \
    --jobs 32 \
    --default-resources mem_mb=8000 runtime=120
```

### Potential future parallelism improvements

The following changes would unlock additional concurrency beyond what Snakemake currently exploits:

- **BLAST per-feature rules** — splitting the Stage 4 BLAST loop into one rule per feature file (see [Optimisation opportunities](#idoptimise)) would add tens of independent jobs to the DAG in parallel with other Stage 4 analysis.
- **Internal R parallelism** — R scripts that iterate over strains or genomic windows could use [`future`](https://future.futureverse.org/) or [`BiocParallel`](https://bioconductor.org/packages/BiocParallel/) to use multiple cores within a single rule, complementing Snakemake-level parallelism.
- **Scatter/gather for Stage 2** — the 13 category fingerprint scripts already run as separate rules per category. Adding a gather rule that merges outputs in parallel (instead of sequentially ordering them) would reduce the critical path through Stage 2.
- **Per-strain report rendering** — the final RMarkdown report currently processes all strains in one R session. Rendering a lightweight per-strain sub-report first (parallelised across strains) and then merging into the final dashboard would reduce the report generation bottleneck.

<br>

[Back to index](#idindex)

<br>

## 8. Visual summary<a name="idsummary"></a>

<img width="1026" height="897" alt="Di-GRAPH" src="https://github.com/acb-lab/Di_GRAPH/blob/b5fc94f42fa9bbaba9324734c433b23ee4b31c6c/images/Di-GRAPH_visual_summary.png" />


[Back to index](#idindex)

