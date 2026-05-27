"""
Pydantic models for Di-GRAPH configuration.

All parameters that were hardcoded in Di-GRAPH_v0.4.0_COMPLETE.sh are now
declared here with types, defaults, and validators. Config is loaded from a
YAML file via ``load_config()``, which resolves relative paths against the
directory that contains the config file.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


# ---------------------------------------------------------------------------
# Sub-models
# ---------------------------------------------------------------------------


class GenomeConfig(BaseModel):
    """Paths to reference genome resources and derived index prefixes."""

    genome_fasta: Path
    """Full PMV reference genome FASTA (e.g. files/RG/RG_PMV_v9.fasta)."""

    genome_fasta_chriii: Path
    """CHRIII-only FASTA used for BWA mutagenic analysis."""

    bowtie_index_prefix: Path
    """Prefix for bowtie/bowtie2 indices; built at runtime if absent."""

    chrom_order: Path
    """Chromosome order file created from the FASTA index."""

    genome_size: int = 14272230
    """Effective genome size (bp) used for RPGC normalisation."""

    @field_validator("genome_fasta", "genome_fasta_chriii", mode="before")
    @classmethod
    def _path_must_exist(cls, v: object) -> Path:
        """Validate that FASTA files exist on disk."""
        p = Path(str(v))
        if not p.exists():
            raise ValueError(f"File not found: {p}")
        return p


class MATCoordinates(BaseModel):
    """Genomic coordinates of the two MAT loci and the HO cut site."""

    chriii_start: int = 199953
    chriii_end: int = 201553
    chrv_start: int = 289025
    chrv_end: int = 290625
    ho_site: int = 200753
    """Position of the HO endonuclease cut on CHRIII."""
    ho_site_upstream: int = 200689
    """Upstream polymorphic position used as a repair-pathway marker."""


class PolymorphismConfig(BaseModel):
    """
    Paired polymorphic positions used for MAT quantification.

    Each ``chriii_positions[i]`` is compared against ``chrv_positions[i]``.
    Coverage at each site is baseline-corrected using the mean of the
    ``baseline_offset`` flanking positions on each side.
    """

    chriii_positions: list[int] = [
        200119, 200167, 200212, 200272, 200326, 200386, 200449,
        200509, 200542, 200575, 200635, 200689, 200753, 200817,
        200882, 200947, 201012, 201077, 201148, 201207, 201272,
        201337, 201402,
    ]
    chrv_positions: list[int] = [
        289191, 289239, 289284, 289344, 289398, 289458, 289521,
        289581, 289614, 289647, 289707, 289761, 289825, 289889,
        289954, 290019, 290084, 290149, 290220, 290278, 290343,
        290408, 290473,
    ]
    baseline_offset: int = 19
    """Lines above/below each target position used to compute baseline coverage."""

    @model_validator(mode="after")
    def _positions_must_be_paired(self) -> "PolymorphismConfig":
        """CHRIII and CHRV coordinate lists must have the same length."""
        if len(self.chriii_positions) != len(self.chrv_positions):
            raise ValueError(
                "chriii_positions and chrv_positions must have equal length "
                f"(got {len(self.chriii_positions)} vs {len(self.chrv_positions)})"
            )
        return self


class TrimConfig(BaseModel):
    """Read trimming and quality-filter parameters."""

    read_length_75: int = 75
    """Length of the 75nt fragments extracted from each paired-end read."""

    read_length_18: int = 18
    """Length of the 18nt fragments used for MAT locus quantification."""

    quality_threshold: Annotated[int, Field(ge=1, le=40)] = 30
    """Minimum per-base Phred quality score (fastp -q flag)."""


class SampleConfig(BaseModel):
    """
    One strain/sample group — a subdirectory inside the working directory
    that contains the paired FASTQ.gz files.
    """

    name: str
    """Subdirectory name, e.g. ``1_Wt``."""

    path: Path
    """Absolute path to the strain directory."""

    is_reference: bool = False
    """Set to ``true`` for the reference strain used in rate comparisons."""

    @field_validator("path", mode="before")
    @classmethod
    def _dir_must_exist(cls, v: object) -> Path:
        p = Path(str(v))
        if not p.is_dir():
            raise ValueError(f"Sample directory not found: {p}")
        return p


class ExperimentConfig(BaseModel):
    """Replicate labels and timepoint names."""

    names: list[str]
    """Experiment/replicate identifiers, e.g. ``[E1, E2, E3]``."""

    timepoints: list[str] = ["T0", "TSG", "TLG", "TLR"]
    """Timepoint prefixes used in FASTQ file names."""

    @field_validator("names")
    @classmethod
    def _at_least_one(cls, v: list[str]) -> list[str]:
        if len(v) < 1:
            raise ValueError("At least one experiment name is required")
        return v


class PathsConfig(BaseModel):
    """All external directory references, resolvable as relative or absolute paths."""

    working_dir: Path
    """Root directory that contains the strain subdirectories."""

    genome_dir: Path
    """Directory with reference genome FASTA files."""

    categories_dir: Path
    """Directory with the 13 genomic category TSV annotation files."""

    blast_dir: Path
    """Directory with per-feature FASTA files used for BLAST cross-validation."""

    report_dir: Path
    """Directory with the RMarkdown template and report generation script."""

    scripts_dir: Path
    """Directory with ``R/`` and ``python/`` sub-folders containing analysis scripts."""

    output_dir: Path | None = None
    """Output directory; defaults to ``working_dir`` if not provided."""

    @model_validator(mode="after")
    def _dirs_must_exist(self) -> "PathsConfig":
        for field in ("working_dir", "genome_dir", "categories_dir", "blast_dir", "report_dir"):
            v: Path = getattr(self, field)
            if not v.is_dir():
                raise ValueError(f"paths.{field}: directory not found: {v}")
        if self.output_dir is None:
            self.output_dir = self.working_dir
        return self


class ResourceConfig(BaseModel):
    """Compute resource settings passed to Snakemake."""

    threads: Annotated[int, Field(ge=1, le=256)] = 4
    """Threads passed to individual tools (bowtie -p, fastp --thread, etc.)."""

    snakemake_cores: Annotated[int, Field(ge=1)] = 4
    """Total cores available to the Snakemake scheduler (--cores)."""

    use_conda: bool = True
    """Whether to activate conda environments per rule (--use-conda)."""


# ---------------------------------------------------------------------------
# Top-level model
# ---------------------------------------------------------------------------


class DiGraphConfig(BaseModel):
    """
    Top-level configuration model for Di-GRAPH.

    Instantiated by ``load_config(path)`` which reads a YAML file,
    resolves relative paths, and validates all fields.
    """

    paths: PathsConfig
    genome: GenomeConfig
    mat_coordinates: MATCoordinates = MATCoordinates()
    polymorphisms: PolymorphismConfig = PolymorphismConfig()
    trimming: TrimConfig = TrimConfig()
    experiments: ExperimentConfig
    samples: list[SampleConfig]
    resources: ResourceConfig = ResourceConfig()

    @model_validator(mode="after")
    def _at_most_one_reference(self) -> "DiGraphConfig":
        """Only one sample can be designated as the reference strain."""
        refs = [s for s in self.samples if s.is_reference]
        if len(refs) > 1:
            raise ValueError(
                f"Only one sample can be the reference strain; got {[s.name for s in refs]}"
            )
        return self

    @property
    def reference_strain(self) -> SampleConfig | None:
        """Return the reference strain if one is designated, else None."""
        for s in self.samples:
            if s.is_reference:
                return s
        return None

    @property
    def sample_names(self) -> list[str]:
        """Convenience accessor for sample name list."""
        return [s.name for s in self.samples]


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


def load_config(path: Path) -> DiGraphConfig:
    """
    Read ``path`` as YAML and return a validated :class:`DiGraphConfig`.

    All relative paths inside the YAML are resolved against the directory
    that contains the config file, making the config portable across machines.

    Args:
        path: Absolute or relative path to ``config.yaml``.

    Returns:
        Fully validated :class:`DiGraphConfig` instance.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        pydantic.ValidationError: If any field fails validation.
    """
    path = Path(path).resolve()
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    base_dir = path.parent
    raw: dict = yaml.safe_load(path.read_text())

    # Resolve all path-like leaf values relative to base_dir before Pydantic
    # validates them.  This allows users to write ``files/RG/...`` instead of
    # absolute paths while still being machine-agnostic.
    raw = _resolve_paths(raw, base_dir)

    return DiGraphConfig.model_validate(raw)


def _resolve_paths(obj: object, base: Path) -> object:
    """
    Recursively walk a YAML-parsed dict and resolve string values that look
    like paths against ``base``.  Only strings containing ``/`` or starting
    with ``.`` are treated as paths; plain strings such as ``E1`` are left
    untouched.
    """
    if isinstance(obj, dict):
        return {k: _resolve_paths(v, base) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_resolve_paths(item, base) for item in obj]
    if isinstance(obj, str) and ("/" in obj or obj.startswith(".")):
        candidate = Path(obj)
        if not candidate.is_absolute():
            return str(base / candidate)
    return obj
