# Configuration

Di-GRAPH is configured through a single YAML file (`config/config.yaml`). All paths can be absolute or relative to the directory that contains the config file, making runs fully portable across machines.

Run `digraph validate --config config/config.yaml` at any time to check the file for errors before launching the pipeline.

---

## `paths`

Directories used by the pipeline. The only value you **must** change is `working_dir`; all others default to paths within the repository.

| Key | Required | Description |
| --- | --- | --- |
| `working_dir` | **yes** | Root directory containing per-strain subdirectories with FASTQ.gz files |
| `genome_dir` | no | Reference genome directory (default: `files/RG/`) |
| `categories_dir` | no | Genomic category TSV annotation files (default: `files/Categories/`) |
| `blast_dir` | no | Per-feature FASTA files for BLAST cross-validation (default: `files/BLAST/features_extraction/`) |
| `report_dir` | no | RMarkdown template directory (default: `files/Report_files/`) |
| `scripts_dir` | no | Analysis scripts (default: `scripts/`) |
| `output_dir` | no | Output root; defaults to `working_dir` if omitted |

```yaml
paths:
  working_dir: /data/experiments/run1
  genome_dir:  files/RG               # relative to this config file
```

---

## `genome`

Reference genome resources for the PMV strain. These values should not need changing unless adapting the pipeline to a different genetic background.

| Key | Default | Description |
| --- | --- | --- |
| `genome_fasta` | `files/RG/RG_PMV_v9.fasta` | Full reference genome FASTA |
| `genome_fasta_chriii` | `files/RG/RG_PMV_v9_CHRIII.fasta` | CHRIII-only FASTA used for BWA mutagenic analysis |
| `bowtie_index_prefix` | `files/RG/S_cerevisiae_indexed` | Prefix for bowtie/bowtie2 indices (built at runtime if absent) |
| `chrom_order` | `files/RG/chrom_order.txt` | Chromosome order file for sorted bedGraph/WIG outputs |
| `genome_size` | `14272230` | Effective genome size in bp used for RPGC normalisation |

---

## `mat_coordinates`

Genomic coordinates of the two MAT loci and the HO cut site. Change only if adapting to a different genetic background.

| Key | Default | Description |
| --- | --- | --- |
| `chriii_start` | `199953` | Start of the *MATa* locus on CHRIII |
| `chriii_end` | `201553` | End of the *MATa* locus on CHRIII |
| `chrv_start` | `289025` | Start of the *MATa′* locus on CHRV |
| `chrv_end` | `290625` | End of the *MATa′* locus on CHRV |
| `ho_site` | `200753` | HO endonuclease cut site on CHRIII |
| `ho_site_upstream` | `200689` | Upstream polymorphic position used as a repair-pathway marker |

---

## `polymorphisms`

Twenty-three paired polymorphic positions used to quantify gene conversion between the *MATa* and *MATa′* loci. Each `chriii_positions[i]` is compared against `chrv_positions[i]`. Coverage at each site is baseline-corrected using the mean of `baseline_offset` flanking positions.

```yaml
polymorphisms:
  chriii_positions: [200119, 200167, ..., 201402]   # 23 values
  chrv_positions:   [289191, 289239, ..., 290473]   # 23 values (paired)
  baseline_offset: 19
```

!!! warning
    `chriii_positions` and `chrv_positions` must have the same length. Validation fails at startup if they differ.

---

## `trimming`

Read trimming and quality-filtering parameters.

| Key | Default | Description |
| --- | --- | --- |
| `read_length_75` | `75` | Length of fragments extracted for coverage and discordant analysis |
| `read_length_18` | `18` | Length of fragments used for MAT locus quantification |
| `quality_threshold` | `30` | Minimum Phred quality score (fastp `-q` flag) |

---

## `experiments`

Replicate and timepoint labels. Values must exactly match the FASTQ file naming convention `<timepoint>_<replicate>_R1.fastq.gz`.

| Key | Default | Description |
| --- | --- | --- |
| `names` | `[E1, E2, E3]` | Replicate identifiers |
| `timepoints` | `[T0, TSG, TLG, TLR]` | Timepoint prefixes |

**Timepoint meanings:**

- **T0** — prior to DSB induction
- **TSG** — Short Galactose (non-selected survivors)
- **TLG** — Long Galactose (selected survivors)
- **TLR** — Long Raffinose (undamaged cells)

---

## `samples`

List of strains to analyse. Each entry must have a `name` matching a subdirectory of `working_dir`.

| Key | Required | Description |
| --- | --- | --- |
| `name` | **yes** | Subdirectory name inside `working_dir` |
| `path` | no | Absolute path to the strain directory (inferred from `working_dir/name` if omitted) |
| `is_reference` | no | Set `true` for the strain used as reference in rate comparisons (at most one) |

```yaml
samples:
  - name: 1_Wt
    is_reference: true
  - name: 2_exo1
  - name: 3_sgs1
```

---

## `resources`

Compute resource settings.

| Key | Default | Description |
| --- | --- | --- |
| `threads` | `4` | Threads passed to individual tools (bowtie `-p`, fastp `--thread`, etc.) |
| `snakemake_cores` | `4` | Total cores available to the Snakemake scheduler (`--cores`) |
| `use_conda` | `true` | Activate conda environments per rule (`--use-conda`) |

!!! tip
    Set `snakemake_cores` to the number of physical cores on your machine. With multiple strains, Snakemake will automatically run up to that many alignment jobs in parallel.
