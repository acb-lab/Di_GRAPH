# Stage 3 — Mutagenic Rate at HO

Stage 3 quantifies the mutagenic signature at the HO cut site and classifies repair pathway choice by examining the nucleotide identity at two diagnostic positions flanking the cut.

---

## Overview

```
.stage1_done
    │
    ├─ bwa mem (T0, TLG, TLR × E1/E2/E3) ──► filter CHRIII ──► extract HO-site reads
    │
    ├─ MUT_extract_200689.py ──► reads_A.bam, reads_G.bam
    │
    ├─ MUT_extract_200753.py (on reads_A) ──► reads_A_T.bam, reads_A_C.bam
    │
    ├─ merge BAM classes ──► bcftools mpileup + igvtools ──► WIG files per variant class
    │
    ├─ MUT_calc_repair_pathway.R ──► repair pathway fractions TSV
    │
    └─ MUT_plot_repair_comparison.R ──► comparison plots ──► .stage3_done
```

---

## Inputs

- `.stage1_done` sentinel.
- CHRIII-only reference genome FASTA (`genome.genome_fasta_chriii`).
- Reads from Stage 1 75 nt filtered FASTQ (T0, TLG, TLR timepoints only).

## Outputs

Per strain:

- `*_HOs.bam` — reads overlapping the HO cut site (CHRIII:200753 ± flanking).
- `*_200689_A.bam`, `*_200689_G.bam` — reads classified by nucleotide at position 200689.
- `*_200689_A_200753_T.bam`, `*_200689_A_200753_C.bam` — reads further classified at 200753.
- WIG coverage files for each of the five variant classes.
- Repair pathway fraction TSVs and comparison SVG plots.
- `.stage3_done` sentinel.

---

## Diagnostic positions

The repair pathway is inferred from the nucleotide observed at two positions near the HO site:

| Position | Nucleotide | Interpretation |
| --- | --- | --- |
| CHRIII:200689 | **A** | wild-type *MATa* sequence (or gene conversion from *MATa′*) |
| CHRIII:200689 | **G** | NHEJ / other mutagenic repair |
| CHRIII:200753 | **T** | *MATa* allele |
| CHRIII:200753 | **C** | *MATa′* allele (gene conversion) |

---

## Variant classes

| Class | Meaning |
| --- | --- |
| `A_CT` | Reads with A at 200689 (any allele at 200753) |
| `A_CT_G` | All reads (A_CT + G at 200689) |
| `A_200753_T` | Gene conversion with *MATa* identity retained |
| `A_200753_C` | Gene conversion incorporating *MATa′* sequence |
| `G` | Reads with G at 200689 (NHEJ / mutagenic events) |

---

## Tools

| Tool / Script | Role |
| --- | --- |
| bwa mem | Alignment of 75 nt reads to the CHRIII reference |
| samtools | BAM filtering, sorting, indexing, merging |
| `MUT_extract_200689.py` | Classify reads by nucleotide at CHRIII:200689 |
| `MUT_extract_200753.py` | Classify reads by nucleotide at CHRIII:200753 |
| bcftools mpileup | Variant-aware pileup for WIG generation |
| igvtools | WIG coverage track generation per variant class |
| `MUT_calc_repair_pathway.R` | Compute repair pathway fractions |
| `MUT_plot_repair_comparison.R` | Multi-strain comparison plots |

!!! note "Timepoints"
    Stage 3 uses only **T0**, **TLG**, and **TLR**. TSG is excluded because the short galactose exposure does not allow sufficient time for repair pathway classification.
