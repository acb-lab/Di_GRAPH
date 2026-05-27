# Stage 2 — Genomic Categories

Stage 2 characterises how coverage changes across 13 genomic feature classes between pairwise timepoint comparisons, providing a genome-wide view of how different genomic elements are affected by DSB induction and repair.

---

## Overview

```
.stage1_done
    │
    ├─ SR_process_categories.R  (× 13 categories, per strain) ──► fingerprint TSVs
    │
    ├─ SR_order_categories.R    (× 3 suffixes: T0vsTSG, T0vsTLG, T0vsTLR) ──► ordered TSVs
    │
    ├─ SR_plot_categories.R     (× 3 suffixes) ──► SVG plots
    │
    └─ SR_plot_gal_vs_raf.R ──► galactose vs raffinose comparison SVG ──► .stage2_done
```

---

## Inputs

- `.stage1_done` sentinel (all Stage 1 rules complete for the strain).
- Coverage TSV produced by Stage 1 (`*_Coverage.tsv`).
- Genomic category annotation TSV files from `paths.categories_dir`.

## Outputs

Per strain:

- `*_<category>_fingerprint.tsv` (× 13) — per-category average coverage fingerprints.
- `*_T0vsTSG_ordered.tsv`, `*_T0vsTLG_ordered.tsv`, `*_T0vsTLR_ordered.tsv` — coverage ordered across timepoints.
- SVG plots for each pairwise comparison and the galactose vs raffinose overview.
- All outputs moved to `Category_Data_TSV_and_Plots/` subdirectory.
- `.stage2_done` sentinel.

---

## Genomic categories

| Category | Description |
| --- | --- |
| ORF | Open reading frames |
| LTR | Long terminal repeats |
| TEG | Ty element genes |
| Ty | Ty retrotransposons |
| tRNA | Transfer RNA genes |
| rRNA | Ribosomal RNA genes |
| ncRNA | Non-coding RNA genes |
| snRNA | Small nuclear RNA genes |
| snoRNA | Small nucleolar RNA genes |
| ARS | Autonomously replicating sequences |
| Cen | Centromeric regions |
| Tel | Subtelomeric regions |
| Int | Intergenic regions |

---

## Pairwise comparisons

| Suffix | Meaning |
| --- | --- |
| `T0vsTSG` | Before DSB induction vs non-selected survivors |
| `T0vsTLG` | Before DSB induction vs selected survivors |
| `T0vsTLR` | Before DSB induction vs undamaged cells |

---

## Tools

| Script | Role |
| --- | --- |
| `SR_process_categories.R` | Compute per-category average coverage fingerprint for one category |
| `SR_order_categories.R` | Order and align fingerprints across timepoints for one comparison suffix |
| `SR_plot_categories.R` | Generate SVG comparison plots for one suffix |
| `SR_plot_gal_vs_raf.R` | Galactose vs raffinose overview comparison |
