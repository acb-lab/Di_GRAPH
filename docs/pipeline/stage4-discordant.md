# Stage 4 — Discordant Reads

Stage 4 identifies and characterises inter-chromosomal discordant read pairs — reads whose two mates map to different chromosomes — as a proxy for gross chromosomal rearrangements. BLAST cross-validation filters spurious discordant signals arising from repetitive elements.

---

## Overview

Stage 4 is divided into four analytical groups:

### Group A — Alignment and discordant extraction

```
75 nt paired FASTQ ──► cutadapt (paired) ──► fastp ──► bowtie (single-read) ──► common mapped reads
                                                     ──► bowtie2 (paired) ──► concordant / discordant BAMs
                                                                           ──► inter-discordant pairs TSV
```

### Group B — MAT coverage from discordant reads

```
inter-discordant pairs ──► extract MAT-overlapping reads ──► bamCoverage ──► DISC_plot_MAT_coverage.R
```

### Group C — 18 nt discordant (TSG and TLG only)

Same trim→filter→align chain applied to 18 nt fragments for 10 read-pair subsets.

### Group D — BLAST cross-validation and R analysis

```
inter-discordant pairs ──► DISC_process_inter_discordant_10kb.R ──► processed pairs TSV
                       ──► blastn (per-feature FASTA, shell loop) ──► BLAST results
                       ──► validation + classification R scripts (≈ 20 rules)
                       ──► hotspots, recombination rates, radar plots, network ──► .stage4_done
```

---

## Inputs

- `.stage1_done` sentinel (reuses 75 nt and 18 nt filtered reads from Stage 1).
- BLAST per-feature FASTA files from `paths.blast_dir`.

## Outputs

Per strain:

- `*_inter_discordant_pairs.tsv` — raw inter-chromosomal discordant pairs.
- `*_inter_discordant_pairs_unique_processed_blast_chromosome.tsv` — BLAST-validated pairs.
- Level-stratified discordant pair tables (option1–option5).
- Coverage tracks at MAT loci derived from discordant reads.
- Global distribution, category distribution, matrix, hotspot, recombination rate, radar, and network SVG plots.
- `.stage4_done` sentinel.

---

## BLAST cross-validation

Each inter-discordant pair is cross-validated against a set of per-feature FASTA databases (Ty elements, LTR sequences, subtelomeric repeats, etc.) to identify pairs that map discordantly due to sequence similarity rather than a true chromosomal rearrangement. Pairs that are explained by BLAST hits are flagged and excluded from the final rearrangement calls.

The BLAST loop runs one `blastn` call per feature file inside the `blast_dir` directory. See the [Optimisation opportunities](../index.md) section for a note on how this could be parallelised in future.

---

## Concordant/discordant classification

Reads are classified by mapping orientation after bowtie2 paired-end alignment. The awk pipeline extracts only pairs where:

- Both mates map with exactly 75 M (no indels or soft-clipping).
- The pair is either concordant (same chromosome, expected orientation) or inter-chromosomal discordant (mates map to different chromosomes).

---

## Tools and scripts

| Tool / Script | Role |
| --- | --- |
| bowtie | Single-read alignment for common-mapped-reads identification |
| bowtie2 | Paired-end alignment for concordant/discordant classification |
| samtools | BAM manipulation and indexing |
| blastn (BLAST 2.16) | Cross-validation of discordant pairs against repeat databases |
| `DISC_process_inter_discordant_10kb.R` | Process and filter discordant pairs within 10 kb windows |
| `DISC_validate_inter_discordant_10kb.R` | Validate BLAST cross-validation results |
| `DISC_global_distribution.R` | Genome-wide distribution of rearrangements |
| `DISC_identify_hotspots.R` | Identify recombination hotspots |
| `DISC_recombination_rate.R` | Calculate per-region recombination rates |
| `DISC_radar_*.R` | Radar chart visualisations |
| `DISC_discordant_network.R` | Network graph of inter-chromosomal interactions |
