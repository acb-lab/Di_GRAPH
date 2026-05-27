# Stage 5 — HTML Report

Stage 5 renders a self-contained HTML flexdashboard report summarising the outputs of all previous stages across every strain.

---

## Overview

```
.stage1_done (all strains)
.stage2_done (all strains)
.stage3_done (all strains)
.stage4_done (all strains)
    │
    └─ generate_report.R ──► Di-GRAPH_report.html
```

The report is generated only after all four analysis stages have completed for every strain in the config.

---

## Inputs

- Stage sentinels `.stage1_done` through `.stage4_done` for all configured strains.
- `Di-GRAPH_report.Rmd` — RMarkdown flexdashboard template (`paths.report_dir`).
- All TSV and SVG outputs produced by Stages 1–4.

## Output

- `<working_dir>/Di-GRAPH_report.html` — standalone HTML report (no server required, opens in any browser).

---

## Report sections

| Section | Contents |
| --- | --- |
| **Overview** | Alignment statistics and run log for each strain |
| **MAT analysis** | Coverage profiles at *MATa* / *MATa′* loci, gene conversion patterns, and polymorphism incorporation |
| **Mutagenic profiling at HO site** | Repair pathway choice fractions and mutagenic signature plots |
| **Genome-wide analysis** | Genomic category fingerprints, discordant read distribution, hotspot maps, and recombination network |

---

## Tools

| Tool | Role |
| --- | --- |
| R 4.4 | Report rendering runtime |
| rmarkdown / knitr | Markdown + R code chunk evaluation |
| flexdashboard 0.6 | HTML dashboard layout |
| `generate_report.R` | Entry-point script: calls `rmarkdown::render()` on the `.Rmd` template |

!!! note
    The report is rendered from the `report_dir` directory so that `Di-GRAPH_report.Rmd` can locate its assets by relative path. The final HTML file is then moved to the `working_dir` root.
