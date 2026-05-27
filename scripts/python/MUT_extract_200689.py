#!/usr/bin/env python3
"""Extract reads with A or G at CHRIII position 200689."""

from pathlib import Path

import pysam
import typer

app = typer.Typer(add_completion=False)


def _extract_reads(bam_file: str, nucleotide: str, output_bam: str, ref_name: str, position: int) -> None:
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
