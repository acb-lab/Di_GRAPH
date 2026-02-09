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

This pipeline is designed to analyze how stage-specific HR repair factors differentially control DSB repair fidelity and genome-wide stability in the *S. cerevisiae* PMV genetic background, which allows the induction of a DSB in the *MATa* *locus* on chromosome III in a galactose-dependent manner; and contains an engineered *MATa'* *locus* on chromosome V used as donor for recombination. 

> Note: detailed information about the PMV genetic background is available in *Ramos et al.,2022 - Cell Reports*, <https://doi.org/10.1016/j.celrep.2021.110201>


By comparing data prior to DSB induction, during its repair (non-selected survivors) and from survivor populations (selected survivors), Di-GRAPH enables the classification of lethal *vs* non-lethal rearrangements arising during the repair of the DSB. Additionally, Di-GRAPH incorporates the assessment of undamaged contitions (undamaged cells) to distinguish DNA damage-dependent from cell-cycle-dependent genomic alterations.

After providing paired-end genomic sequencing data from these 4 different timepoints, Di-GRAPH will:

- i)  Perform *MATa/MATa'* *loci* coverage profiling to identify gene conversion patterns and polymporphisms incorporation.

- ii) Characterize the frequency, directionality and extent of individual gene conversion events between *MATa/MATa'* *loci* by applying inter-chromosomal discordant read mapping.

- iii) Assess HO associated mutagenic pattern and repair pathway choice upon DSB induction in the *MATa* *locus*.

- iv) Evaluate global genome stability by characterizing how different genomic categories (e.g. ORFs, intergenic regions, LTR, TEG, Ty, tRNA, rRNA, ncRNA, snRNA, snoRNA, ARS, centromere and subtelomeric regions) behave in cells lacking defined HR factors, both after DSB induction and in the absence of DNA damage.

- v) Apply inter-chromosomal discordant read mapping combined with BLAST cross-validation to identify and map genome-wide chromosomal rearrangements in both DNA-damaged and undamaged cells, characterizing how the genome is reshaped in the absence of defined HR factors.

<br>

You can find a **visual summary** of the different steps conducted by Di-GRAPH [here](#idsummary).

<br>

[Back to index](#idindex)

<br>

## 2. Installation<a name="idinstall"></a>

Di-GRAPH is intended to be run by using conda, which you can install by following [this instructions](https://www.anaconda.com/docs/getting-started/miniconda/install).  

Once you have conda installed, you can download Di-GRAPH repository:

```
conda install git
git clone https://github.com/acb-lab/Di_GRAPH.git
```

An install a conda environment (`digraph`) with all the required software and dependencies:

```
conda env create -f Di_GRAPH/environment/digraph.yml -n digraph
```

Activate `digraph` environment and install Di-GRAPH executable 

```
conda activate digraph

# Make Di-GRAPH.sh executable
chmod 700 Di_GRAPH/Di-GRAPH.sh # 770 for group permissions, 777 for all permissions

# Add Di-GRAPH.sh and required files to PATH
cp Di_GRAPH/Di-GRAPH.sh $CONDA_PREFIX/bin/Di-GRAPH.sh
cp -r Di_GRAPH/files $CONDA_PREFIX/files
```

\* *We are working to integrate everything as a single conda package*

<br>

[Back to index](#idindex)

<br>

## 3. Instructions<a name="idinstr"></a>

Di-GRAPH is designed to be run as a single bash script `Di-GRAPH.sh`. As noted in the [Introduction](#idintro), Di-GRAPH requires paired-end genomic data (FASTQ files) from 4 different timepoints to work. To ensure HO associated mutagenic pattern is performed, paired-end reads should be of a minimum length of 150bp. In addition, Di-GRAPH requires three independent replicates of each timepoint to perform the analysis.

Di-GRAPH also requires the PMV reference genome, the PMV genomic categories and features annotation files, a fasta genomic database for BLAST cross-validation and a folder containing the template for the HTML output report to be used for the analysis. All these files are available in the `Di_GRAPH/files` folder. 

Di-GRAPH will require the user to define the path to the mentioned files with the following options:

```
OBLIGATORY OPTIONS:
    -b/--blast          Path to BLAST genomic database for cross-validation directory (provided in GitHub)
    -c/--categories     Path to genomic categories and features annotation directory (provided in GitHub)
    -g/--genome         Path to reference genome directory (provided in GitHub)
    -i/--input          Path to input directory (working directory generated by the user)
    -r/--report         Path to report directory (provided in GitHub)

OTHER OPTIONS:
    -t/--threads        Number of threads to use (default=1)
    -v/--version        Show version
```

<br> 


Di-GRAPH requires the following folder/subfolder infrastructure to perform the complete analysis:

    ------- Folder generated by the user -------
    Working directory folder (working_directory)
    ├── Strain_1 subfolder
    │   ├── T0_E1_R1.fastq.gz
    │   ├── T0_E1_R2.fastq.gz
    │   ├── TSG_E1_R1.fastq.gz
    │   ├── TSG_E1_R2.fastq.gz
    │   ├── ...
    ├── Strain_2 subfolder
    │   ├── .fastq.gz files
    └── Strain_n subfolder
        ├── .fastq.gz files

    ------- Folders available in the Di_GRAPH/files folder from this repository -------
    Reference genome folder (RG)
    ├── RG_PMV_v9.fasta
    ├── ...

    Genomic categories and features annotation folder (Categories)
    ├── PMV_categories.tsv
    ├── 1.PMV.ORF.tsv     
    └── ...

    BLAST folder (BLAST)
    ├── features_extraction
        ├── YLL039C_seq.fasta
        ├── YORWTy1_2_seq.fasta
        ├── ...

    Report folder (Report_files)
    ├── Di-GRAPH_report.Rmd
    ├── ...
 

For each strain to be analyzed, a subfolder inside the working directory folder should be created. Paired-end genomic data (.fastq.gz files) should be placed inside this subfolder and named as follows: `timepoint_replicate_R1.fastq.gz` and `timepoint_replicate_R2.fastq.gz`.

> Note: timepoint names should be defined as **T0** (Timepoint 0, data prior to DSB induction), **TSG** (Timepoint Short Galactose, data from non-selected survivors), **TLG** (Timepoint Long Galactose, data from selected survivors) and **TLR** (Timepoint Long Raffinose, data from undamaged cells). Replicate names should be defined as **E1** (Experiment 1), **E2** (Experiment 2) and **E3** (Experiment 3).
  
<br>

You can see an example of Di-GRAPH usage in [this section](#idexample)

<br>

[Back to index](#idindex)

<br>

## 4. Example usage<a name="idexample"></a>

You can test Di-GRAPH with paired-end genomic data from wild-type PMV cells available in the `Di-GRAPH wild-type dataset` Zenodo repository (https://doi.org/10.5281/zenodo.18470248).  

```bash
## Generate working directory
mkdir working_directory
cd working_directory

## Download wild-type dataset from Zenodo
wget https://zenodo.org/records/18470248/files/1_Wt.tgz
tar xvzf 1_Wt.tgz
rm 1_Wt.tgz

```

Note that for running the current version of Di-GRAPH (we are working to ease this process), you need to maintain the files infrastructure as in [Instructions](#idinstr).  

Then, you can run Di-GRAPH by simply typing:

```bash
conda activate digraph
Di-GRAPH.sh -b $CONDA_PREFIX/files/BLAST -c $CONDA_PREFIX/files/Categories -g $CONDA_PREFIX/files/RG -i working_directory/ -r $CONDA_PREFIX/files/Report_files
```

<br>


[Back to index](#idindex)

<br>

## 5. Expected output<a name="idoutput"></a>

For each strain defined in the working directory, Di-GRAPH will perform `bowtie/bowtie2/bwa` genomic alignments, characterize gene conversion products between *MATa/MATa'* *loci*, define HO associated mutagenic pattern, analyze coverage data regarding all genomic categories and identify global genomic rearrangements. The output files and plots will be stored in the `working_directory/Strain_n` subfolder(s).

To facilitate the interpretation of the results, Di-GRAPH will generate a `Di-GRAPH_report.html` file in the `working_directory` folder. This report is a summary of the results obtained by Di-GRAPH for each strain and includes the following sections:
- **Overview:** contains the alignment statistics and the script log file.
- **MAT analysis:** contains coverage analysis and gene conversion analysis for the *MATa/MATa'* *loci*. It also contains the polymorphisms incorporation analysis.
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

