#!/usr/bin/env python3
"""Extract reads with A or G at CHRIII position 200689."""

from pathlib import Path

import pysam
import typer

app = typer.Typer(add_completion=False)


def _extract_reads(bam_file: str, nucleotide: str, output_bam: str, ref_name: str, position: int) -> None:
    """
    Write to ``output_bam`` all reads in ``bam_file`` that carry ``nucleotide``
    at reference position ``position`` on chromosome ``ref_name``.

    The function uses pysam's ``get_aligned_pairs(matches_only=True)`` to map
    each query base back to its reference coordinate; reads that do not overlap
    the target position (e.g. soft-clipped or deleted at that site) are silently
    skipped.

    Args:
        bam_file:   Path to the sorted, indexed input BAM.
        nucleotide: Single character (``"A"``, ``"G"``, ``"T"``, or ``"C"``)
                    to match against the query base.
        output_bam: Path to write the filtered output BAM (header copied from input).
        ref_name:   Reference sequence name to fetch (e.g. ``"CHRIII"``).
        position:   1-based reference position to inspect.

    Raises:
        SystemExit: On any pysam error, prints to stderr and exits with code 1.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        out_bam = pysam.AlignmentFile(output_bam, "wb", template=bam)

        for read in bam.fetch(ref_name, position - 1, position):
            if not read.query_sequence:
                continue
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                if ref_pos == position - 1:
                    if read.query_sequence[query_pos] == nucleotide:
                        out_bam.write(read)
                    break

        bam.close()
        out_bam.close()

    except Exception as e:
        typer.echo(f"ERROR processing {bam_file} ({nucleotide}): {e}", err=True)
        raise typer.Exit(code=1)


@app.command()
def main(
    input_bam: Path = typer.Argument(..., help="Input BAM file."),
    output_a: Path = typer.Argument(..., help="Output BAM — reads with A at CHRIII:200689."),
    output_g: Path = typer.Argument(..., help="Output BAM — reads with G at CHRIII:200689."),
) -> None:
    """Classify HO-site reads by nucleotide at CHRIII:200689 (A or G)."""
    _extract_reads(str(input_bam), "A", str(output_a), "CHRIII", 200689)
    _extract_reads(str(input_bam), "G", str(output_g), "CHRIII", 200689)


if __name__ == "__main__":
    app()
