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
6.  [Visual summary](#idsummary)

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

## 6. Visual summary<a name="idsummary"></a>

<img width="1026" height="897" alt="Di-GRAPH" src="https://github.com/acb-lab/Di_GRAPH/blob/b5fc94f42fa9bbaba9324734c433b23ee4b31c6c/images/Di-GRAPH_visual_summary.png" />


[Back to index](#idindex)

