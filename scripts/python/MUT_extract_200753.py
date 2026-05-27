#!/usr/bin/env python3
"""Extract reads with T or C at CHRIII position 200753."""

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
    output_t: Path = typer.Argument(..., help="Output BAM — reads with T at CHRIII:200753."),
    output_c: Path = typer.Argument(..., help="Output BAM — reads with C at CHRIII:200753."),
) -> None:
    """Classify HO-site reads by nucleotide at CHRIII:200753 (T or C)."""
    _extract_reads(str(input_bam), "T", str(output_t), "CHRIII", 200753)
    _extract_reads(str(input_bam), "C", str(output_c), "CHRIII", 200753)


if __name__ == "__main__":
    app()
