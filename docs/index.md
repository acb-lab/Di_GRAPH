---
template: home.html
---

## Five-stage pipeline

<div class="dg-stages">
  <a class="dg-stage" href="pipeline/stage1-coverage/">
    <div class="dg-stage__num">01</div>
    <div class="dg-stage__name">Coverage</div>
    <div class="dg-stage__desc">Read trimming, alignment, RPGC-normalised tracks, and MAT locus polymorphism quantification.</div>
  </a>
  <a class="dg-stage" href="pipeline/stage2-categories/">
    <div class="dg-stage__num">02</div>
    <div class="dg-stage__name">Categories</div>
    <div class="dg-stage__desc">Coverage fingerprints across 13 genomic feature classes — ORF, LTR, Ty, tRNA, centromere, and more.</div>
  </a>
  <a class="dg-stage" href="pipeline/stage3-mutagenic/">
    <div class="dg-stage__num">03</div>
    <div class="dg-stage__name">Mutagenic</div>
    <div class="dg-stage__desc">HO-site read extraction, nucleotide classification, and repair pathway choice analysis.</div>
  </a>
  <a class="dg-stage" href="pipeline/stage4-discordant/">
    <div class="dg-stage__num">04</div>
    <div class="dg-stage__name">Discordant</div>
    <div class="dg-stage__desc">Inter-chromosomal discordant mapping, BLAST cross-validation, and rearrangement profiling.</div>
  </a>
  <a class="dg-stage" href="pipeline/stage5-report/">
    <div class="dg-stage__num">05</div>
    <div class="dg-stage__name">Report</div>
    <div class="dg-stage__desc">Self-contained flexdashboard HTML report summarising all stages across every strain.</div>
  </a>
</div>

---

## Quick install

<div class="dg-quickstart">
<div class="dg-quickstart__title">Three commands to get started</div>

```bash
conda env create -f environment/digraph.yml -n digraph && conda activate digraph
pip install -e .
digraph run --config config/config.yaml --cores 8
```

</div>

See the full [installation guide](installation.md) and [configuration reference](configuration.md).

---

## What Di-GRAPH measures

Di-GRAPH is designed for the **PMV genetic background**, which allows galactose-inducible induction
of a single DSB at *MATa* (CHRIII:200753) and provides an engineered *MATa′* donor locus on chromosome V.
By comparing data from four timepoints — prior to break induction, during non-selected and selected
repair, and in undamaged cells — the pipeline classifies:

- Gene conversion frequency, directionality, and extent between *MATa / MATa′* loci.
- Repair pathway choice (gene conversion *vs* NHEJ *vs* mutagenic events) at the HO site.
- Genome-wide genomic stability across 13 feature categories.
- Inter-chromosomal rearrangements and recombination hotspots.

> **Reference:** Ramos *et al.*, 2022 — *Cell Reports* · [doi:10.1016/j.celrep.2021.110201](https://doi.org/10.1016/j.celrep.2021.110201)

---

## Parallelism out of the box

Snakemake automatically runs all independent jobs concurrently.
With 5 strains × 4 timepoints × 3 replicates and `--cores 16`,
**Stage 1 alone dispatches up to 60 alignment jobs in parallel** — no configuration needed.

```bash
digraph run --config config/config.yaml --cores 16   # local multi-core
```

See [Parallelisation](https://github.com/acb-lab/Di_GRAPH#7-parallelisation) for cluster and cloud execution options.
