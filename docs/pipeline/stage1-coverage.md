# Stage 1 — Coverage

Stage 1 performs read trimming, alignment, and coverage track generation for all strains, timepoints, and replicates. It produces the input data consumed by every downstream stage.

---

## Overview

```
FASTQ.gz (R1 + R2)
    │
    ├─ cutadapt (8 × 75 nt windows) ──► filter (fastp) ──► bowtie ──► bamCoverage (RPGC)
    │                                                               ──► Coverage bedGraph
    │
    └─ cutadapt (18 nt) ─────────────► filter (fastp) ──► bowtie ──► Coverage bedGraph (CHRIII + CHRV)
                                                                   ──► Polymorphism diff (Python)
                                                                   ──► Non-RPGC coverage (TSG/TLG/TLR)
    │
    └─► R scripts: SR_process_cov_iii_v.R, SR_process_cov_18nt.R ──► .stage1_done
```

---

## Inputs

- Paired-end FASTQ.gz files named `<timepoint>_<replicate>_R1.fastq.gz` / `R2.fastq.gz` inside each strain directory.
- PMV reference genome FASTA (`genome.genome_fasta`).

## Outputs

Per strain, per timepoint, per replicate:

- `*_75nt_sorted.bam` — sorted, indexed 75 nt alignment BAM.
- `*_RPGC.bedGraph` — RPGC-normalised whole-genome coverage.
- `*_CHRIII_18nt.bedGraph`, `*_CHRV_18nt.bedGraph` — 18 nt coverage at the MAT loci.
- `*_polymorphism_diff.tsv` — baseline-subtracted coverage at each of the 23 polymorphic positions.
- `*_Coverage.tsv` — merged coverage table used by Stage 2.
- `.stage1_done` — sentinel file consumed by Stage 2.

---

## Tools and parameters

| Tool | Version | Role |
| --- | --- | --- |
| cutadapt | 5.1 | Extract 75 nt and 18 nt fragments from 150 bp reads |
| fastp | 0.25 | Quality filter (Phred ≥ `trimming.quality_threshold`) |
| bowtie | 1.3.1 | Short-read alignment (75 nt and 18 nt) |
| samtools | 1.22 | SAM→BAM conversion, sorting, indexing |
| deeptools bamCoverage | 3.5 | RPGC-normalised coverage tracks |
| bedtools | 2.31 | Per-base coverage from sorted BAM |

---

## Key design notes

**75 nt trimming chain** — Eight 75 bp windows are extracted from each 150 bp paired-end read using separate cutadapt calls, then concatenated before alignment. This maximises the number of uniquely-mappable reads while keeping the fragment length consistent for RPGC normalisation.

**18 nt fragments** — 18 nt reads are used exclusively for MAT locus quantification (CHRIII and CHRV) at timepoints TSG, TLG, and TLR. They are not produced for T0, where DSB induction has not occurred.

**Polymorphism diff** — A Python `run:` block in Snakemake reads coverage values at each of the 23 polymorphic positions and subtracts a local baseline (mean of ± `polymorphisms.baseline_offset` flanking positions). This replaces the original awk loop and reads all coordinates from `config.yaml`.

**Bowtie index** — Built once per pipeline run from the FASTA if the index files are absent; subsequent rules reuse the cached index.
