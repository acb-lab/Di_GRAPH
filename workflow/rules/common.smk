"""
Common helpers shared across all stage rule files.

Defines Python convenience functions and wildcard constraints so every rule
file can reference consistent values without duplicating config lookups.
"""

# ---------------------------------------------------------------------------
# Convenience accessors for config values
# ---------------------------------------------------------------------------

def get_sample_names():
    """Return the list of strain/sample names defined in the config."""
    return [s["name"] for s in config["samples"]]


def get_timepoints():
    """Return the list of timepoint labels (T0, TSG, TLG, TLR)."""
    return config["experiments"]["timepoints"]


def get_experiments():
    """Return the list of replicate identifiers (E1, E2, E3 …)."""
    return config["experiments"]["names"]


def get_working_dir():
    """Return the experiment working directory as a string."""
    return config["paths"]["working_dir"]


def get_scripts_dir():
    """Return the scripts directory (contains R/ and python/ sub-dirs)."""
    return config["paths"]["scripts_dir"]


def get_genome_fasta():
    return config["genome"]["genome_fasta"]


def get_genome_fasta_chriii():
    return config["genome"]["genome_fasta_chriii"]


def get_bowtie_index():
    return config["genome"]["bowtie_index_prefix"]


def get_chrom_order():
    return config["genome"]["chrom_order"]


def get_genome_size():
    return config["genome"]["genome_size"]


def get_threads():
    return config["resources"]["threads"]


def get_categories_dir():
    return config["paths"]["categories_dir"]


def get_blast_dir():
    return config["paths"]["blast_dir"]


def get_report_dir():
    return config["paths"]["report_dir"]


# Timepoints that have non-RPGC 18nt coverage (for MAT quantification)
QUANT_TIMEPOINTS = ["TSG", "TLG", "TLR"]

# Timepoints used for mutagenic analysis (BWA / Python scripts)
MUT_TIMEPOINTS = ["T0", "TLG", "TLR"]

# 13 genomic categories
CATEGORIES = ["ORF", "LTR", "TEG", "Ty", "tRNA", "rRNA", "ncRNA",
              "snRNA", "snoRNA", "ARS", "Cen", "Tel", "Int"]

# Polymorphic variant classes used for BCFtools / igvtools
POLY_VARIANTS = ["A_CT", "A_CT_G", "A_200753_T", "A_200753_C", "G"]


# ---------------------------------------------------------------------------
# Wildcard constraints
# Keep wildcards from greedily matching path separators or unrelated strings.
# ---------------------------------------------------------------------------

wildcard_constraints:
    strain     = "[^/]+",
    timepoint  = "T0|TSG|TLG|TLR",
    experiment = "[A-Za-z0-9]+",
    category   = "|".join(CATEGORIES),
    poly       = "|".join(POLY_VARIANTS),
    chr        = "CHRIII|CHRV",
