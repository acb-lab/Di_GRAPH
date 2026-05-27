#!/usr/bin/env python3

#### Python script to extract reads with C or T at postion 200753
### Loop for each strain/sample/experiment in MYWD
### 11/04/2026 - Lydia

import pysam
import sys

def extract_reads(bam_file, nucleotide, output_bam, ref_name, position):
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        out_bam = pysam.AlignmentFile(output_bam, "wb", template=bam)

        for read in bam.fetch(ref_name, position-1, position):
            if not read.query_sequence:
                continue

            aligned_pairs = read.get_aligned_pairs(matches_only=True)

            for query_pos, ref_pos in aligned_pairs:
                if ref_pos == position - 1:
                    if read.query_sequence[query_pos] == nucleotide:
                        out_bam.write(read)
                        break

        bam.close()
        out_bam.close()

    except Exception as e:
        print(f"ERROR processing {bam_file} ({nucleotide}): {e}", file=sys.stderr)
        sys.exit(1)


# ===== Main =====
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: extract_reads.py <input.bam> <output_T.bam> <output_C.bam>", file=sys.stderr)
        sys.exit(1)

    bam_file = sys.argv[1]
    output_bam_T = sys.argv[2]
    output_bam_C = sys.argv[3]

    ref_name = "CHRIII"
    position = 200753

    extract_reads(bam_file, "T", output_bam_T, ref_name, position)
    extract_reads(bam_file, "C", output_bam_C, ref_name, position)