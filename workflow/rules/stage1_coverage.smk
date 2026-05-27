"""
Stage 1 — Coverage and Polymorphisms Analysis.

Workflow:
  1. Build bowtie genome index (once, shared across all samples).
  2. Per strain × timepoint × experiment:
       decompress → trim (75nt + 18nt) → quality filter (Q30)
       → align (bowtie -m1/-k1) → BAM → BedGraph → TSV
       → awk coordinate enrichment → ordered 18nt TSVs
  3. Polymorphism difference calculation (Python run: block).
  4. R scripts for MAT coverage visualisation.

All intermediate files that are not needed downstream are deleted inside the
shell blocks to avoid excessive disk use (mirrors the original bash behaviour).
"""

# ---------------------------------------------------------------------------
# Build bowtie genome index (done once; outputs are .ebwt files)
# ---------------------------------------------------------------------------

rule build_bowtie_index:
    """Build the bowtie1 genome index from the reference FASTA."""
    input:
        fasta = get_genome_fasta(),
    output:
        # bowtie-build produces six .ebwt files; track the first one as a proxy
        idx = get_bowtie_index() + ".1.ebwt",
    params:
        prefix = get_bowtie_index(),
    threads: get_threads()
    log:
        f"{get_working_dir()}/logs/build_bowtie_index.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        "bowtie-build {input.fasta} {params.prefix} > {log} 2>&1"


rule build_fai_and_chrom_order:
    """
    Index the reference FASTA with samtools faidx and extract the chromosome
    order file used by bedtools sort.
    """
    input:
        fasta = get_genome_fasta(),
    output:
        chrom_order = get_chrom_order(),
    log:
        f"{get_working_dir()}/logs/build_fai_chrom_order.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        samtools faidx {input.fasta} 2> {log}
        cut -f1,2 {input.fasta}.fai > {output.chrom_order}
        """


# ---------------------------------------------------------------------------
# Stage 1: per strain × timepoint × experiment
# ---------------------------------------------------------------------------

rule decompress_and_trim_75nt:
    """
    Decompress paired FASTQ.gz, concatenate R1+R2, and produce two 75nt
    fragments per read (head and tail trim via cutadapt).

    Outputs:
      _1.fastq  — last 75 nt of each read  (--cut -75)
      _2.fastq  — first 75 nt of each read (--cut 75)
    """
    input:
        r1 = "{wd}/{strain}/{timepoint}_{experiment}_R1.fastq.gz",
        r2 = "{wd}/{strain}/{timepoint}_{experiment}_R2.fastq.gz",
    output:
        frag1 = temp("{wd}/{strain}/{timepoint}_{experiment}_1.fastq"),
        frag2 = temp("{wd}/{strain}/{timepoint}_{experiment}_2.fastq"),
    params:
        threads = get_threads(),
    log:
        "{wd}/{strain}/logs/stage1_trim75_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        # Decompress (keep originals with -k)
        gunzip -k {input.r1} {input.r2}
        r1={wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_R1.fastq
        r2={wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_R2.fastq

        # Concatenate R1 and R2
        cat "$r1" "$r2" > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}.fastq
        rm "$r1" "$r2"

        # Trim to 75 nt fragments
        cutadapt -j {params.threads} --cut -75 \
            -o {output.frag1} \
            {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}.fastq \
            >> {log} 2>&1
        cutadapt -j {params.threads} --cut 75 \
            -o {output.frag2} \
            {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}.fastq \
            >> {log} 2>&1

        rm {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}.fastq
        """


rule trim_18nt:
    """
    Further trim each 75nt fragment into an 18nt core by two sequential
    cutadapt passes (--cut -38/38, then --cut -19/19), then concatenate
    all eight 18nt outputs into a single file.
    """
    input:
        frag1 = "{wd}/{strain}/{timepoint}_{experiment}_1.fastq",
        frag2 = "{wd}/{strain}/{timepoint}_{experiment}_2.fastq",
    output:
        merged_18nt = temp("{wd}/{strain}/{timepoint}_{experiment}_18nt_nonfiltered.fastq"),
    params:
        threads = get_threads(),
        prefix  = "{wd}/{strain}/{timepoint}_{experiment}",
    log:
        "{wd}/{strain}/logs/stage1_trim18_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        p={params.prefix}

        # Second trim: extract inner 38nt from each 75nt fragment
        cutadapt -j {params.threads} --cut -38 -o ${{p}}_1_1.fastq {input.frag1} >> {log} 2>&1
        cutadapt -j {params.threads} --cut  38 -o ${{p}}_1_2.fastq {input.frag1} >> {log} 2>&1
        cutadapt -j {params.threads} --cut -38 -o ${{p}}_2_1.fastq {input.frag2} >> {log} 2>&1
        cutadapt -j {params.threads} --cut  38 -o ${{p}}_2_2.fastq {input.frag2} >> {log} 2>&1

        # Third trim: extract inner 18nt from each 38nt fragment
        for sub in 1_1 1_2 2_1 2_2; do
            cutadapt -j {params.threads} --cut -19 -o ${{p}}_18nts_${{sub}}_1.fastq ${{p}}_${{sub}}.fastq >> {log} 2>&1
            cutadapt -j {params.threads} --cut  19 -o ${{p}}_18nts_${{sub}}_2.fastq ${{p}}_${{sub}}.fastq >> {log} 2>&1
            rm ${{p}}_${{sub}}.fastq
        done

        # Merge all eight 18nt fragments
        cat ${{p}}_18nts_*.fastq > {output.merged_18nt}
        rm ${{p}}_18nts_*.fastq
        """


rule filter_75nt:
    """Quality-filter the 75nt reads with fastp (Q30)."""
    input:
        frag1 = "{wd}/{strain}/{timepoint}_{experiment}_1.fastq",
        frag2 = "{wd}/{strain}/{timepoint}_{experiment}_2.fastq",
    output:
        filtered = temp("{wd}/{strain}/{timepoint}_{experiment}_75nt.fastq"),
        html     = "{wd}/{strain}/{timepoint}_{experiment}_75nt.fastq.html",
        json     = "{wd}/{strain}/{timepoint}_{experiment}_75nt.fastq.json",
    params:
        threads   = get_threads(),
        q_thresh  = config["trimming"]["quality_threshold"],
    log:
        "{wd}/{strain}/logs/stage1_filter75_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        cat {input.frag1} {input.frag2} > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_75nt_nonfiltered.fastq
        fastp \
            -i {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_75nt_nonfiltered.fastq \
            -o {output.filtered} \
            -q {params.q_thresh} -u 0 -e 30 \
            --thread {params.threads} \
            --html {output.html} \
            --json {output.json} \
            >> {log} 2>&1
        rm {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_75nt_nonfiltered.fastq \
           {input.frag1} {input.frag2}
        """


rule filter_18nt:
    """Quality-filter the 18nt reads with fastp (Q30)."""
    input:
        raw = "{wd}/{strain}/{timepoint}_{experiment}_18nt_nonfiltered.fastq",
    output:
        filtered = temp("{wd}/{strain}/{timepoint}_{experiment}_18nt.fastq"),
        html     = "{wd}/{strain}/{timepoint}_{experiment}_18nt.fastq.html",
        json     = "{wd}/{strain}/{timepoint}_{experiment}_18nt.fastq.json",
    params:
        threads  = get_threads(),
        q_thresh = config["trimming"]["quality_threshold"],
    log:
        "{wd}/{strain}/logs/stage1_filter18_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        fastp \
            -i {input.raw} \
            -o {output.filtered} \
            -q {params.q_thresh} -u 0 -e 30 \
            --thread {params.threads} \
            --html {output.html} \
            --json {output.json} \
            >> {log} 2>&1
        rm {input.raw}
        """


rule align_75nt:
    """
    Align 75nt reads with bowtie1 in unique-match mode (-m 1 -v 0).
    Produces a SAM that is immediately converted to sorted BAM.
    """
    input:
        fastq = "{wd}/{strain}/{timepoint}_{experiment}_75nt.fastq",
        idx   = get_bowtie_index() + ".1.ebwt",
    output:
        bam = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_75nt.bam",
        bai = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_75nt.bai",
    params:
        threads = get_threads(),
        prefix  = get_bowtie_index(),
        sam     = temp("{wd}/{strain}/{timepoint}_{experiment}_75nt.sam"),
    log:
        "{wd}/{strain}/logs/stage1_align75_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        bowtie -p {params.threads} -m 1 -v 0 -S \
            -x {params.prefix} \
            {input.fastq} > {params.sam} 2>> {log}
        samtools sort -@ {params.threads} -o {output.bam} {params.sam}
        samtools index {output.bam} {output.bai}
        rm {params.sam} {input.fastq}
        """


rule align_18nt:
    """
    Align 18nt reads with bowtie1 in multi-match mode (-k 1 -v 0).
    Produces a sorted BAM at the MAT locus region.
    """
    input:
        fastq = "{wd}/{strain}/{timepoint}_{experiment}_18nt.fastq",
        idx   = get_bowtie_index() + ".1.ebwt",
    output:
        bam = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_18nt.bam",
        bai = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_18nt.bai",
    params:
        threads = get_threads(),
        prefix  = get_bowtie_index(),
        sam     = temp("{wd}/{strain}/{timepoint}_{experiment}_18nt.sam"),
    log:
        "{wd}/{strain}/logs/stage1_align18_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        bowtie -p {params.threads} -k 1 -v 0 -S \
            -x {params.prefix} \
            {input.fastq} > {params.sam} 2>> {log}
        samtools sort -@ {params.threads} -o {output.bam} {params.sam}
        samtools index {output.bam} {output.bai}
        rm {params.sam} {input.fastq}
        """


rule coverage_75nt:
    """
    Generate per-base RPGC-normalised BedGraph coverage from 75nt BAM,
    sort it with bedtools, then produce a TSV via bedtools coverage.
    Truncates coverage decimals to 5 places to match original behaviour.
    """
    input:
        bam         = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_75nt.bam",
        bai         = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_75nt.bai",
        chrom_order = get_chrom_order(),
    output:
        tsv = "{wd}/{strain}/Coverage_data/{timepoint}_{experiment}_75nt.tsv",
    params:
        threads     = get_threads(),
        genome_size = get_genome_size(),
        bedgraph    = temp("{wd}/{strain}/{timepoint}_{experiment}_75nt.bedgraph"),
        sorted_bg   = temp("{wd}/{strain}/{timepoint}_{experiment}_75nt_sorted.bedgraph"),
    log:
        "{wd}/{strain}/logs/stage1_cov75_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        bamCoverage -b {input.bam} -o {params.bedgraph} \
            -of bedgraph -p {params.threads} -bs 1 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize {params.genome_size} \
            >> {log} 2>&1

        bedtools sort -i {params.bedgraph} -g {input.chrom_order} > {params.sorted_bg}

        bedtools coverage \
            -a {params.sorted_bg} -b {input.bam} \
            -sorted -g {input.chrom_order} -d \
            > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_75nt_raw.tsv 2>> {log}

        # Truncate coverage column (col 4) to 5 decimal places
        awk -F'\\t' 'BEGIN{{OFS=FS}} {{
            split($4,a,".");
            if(length(a)>1) $4=a[1]"."substr(a[2],1,5);
            print
        }}' {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_75nt_raw.tsv \
        > {output.tsv}

        rm {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_75nt_raw.tsv \
           {params.bedgraph} {params.sorted_bg}
        """


rule extract_chr_coverage_75nt:
    """
    Extract CHRIII and CHRV columns from the 75nt TSV, add absolute
    coordinates (NR), and produce MAT-region subsets.
    """
    input:
        tsv = "{wd}/{strain}/Coverage_data/{timepoint}_{experiment}_75nt.tsv",
    output:
        chriii     = "{wd}/{strain}/{timepoint}_{experiment}_75nt_CHRIII.tsv",
        chrv       = "{wd}/{strain}/{timepoint}_{experiment}_75nt_CHRV.tsv",
        chriii_mat = "{wd}/{strain}/{timepoint}_{experiment}_75nt_CHRIII_MATa.tsv",
        chrv_mat   = "{wd}/{strain}/{timepoint}_{experiment}_75nt_CHRV_MATa.tsv",
    params:
        chriii_start = config["mat_coordinates"]["chriii_start"],
        chriii_end   = config["mat_coordinates"]["chriii_end"],
        chrv_start   = config["mat_coordinates"]["chrv_start"],
        chrv_end     = config["mat_coordinates"]["chrv_end"],
    log:
        "{wd}/{strain}/logs/stage1_extchr_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        awk -F'\\t' '$1=="CHRIII" {{print $1, $4}}' OFS='\\t' {input.tsv} | \
            awk 'BEGIN{{OFS="\\t"}} {{print $0, NR}}' > {output.chriii}

        awk -F'\\t' '$1=="CHRV" {{print $1, $4}}' OFS='\\t' {input.tsv} | \
            awk 'BEGIN{{OFS="\\t"}} {{print $0, NR}}' > {output.chrv}

        awk -F'\\t' '$3 >= {params.chriii_start} && $3 <= {params.chriii_end}' \
            {output.chriii} > {output.chriii_mat}

        awk -F'\\t' '$3 >= {params.chrv_start} && $3 <= {params.chrv_end}' \
            {output.chrv} > {output.chrv_mat}
        """


rule coverage_18nt:
    """
    Generate per-base RPGC-normalised 18nt BedGraph coverage at the two
    MAT loci (CHRIII and CHRV), sort, and produce ordered TSVs with
    absolute genomic coordinates.
    """
    input:
        bam         = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_18nt.bam",
        bai         = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_18nt.bai",
        chrom_order = get_chrom_order(),
    output:
        chriii_tsv  = "{wd}/{strain}/{timepoint}_{experiment}_CHRIII_18nt_ordered.tsv",
        chrv_tsv    = "{wd}/{strain}/{timepoint}_{experiment}_CHRV_18nt_ordered.tsv",
    params:
        threads      = get_threads(),
        genome_size  = get_genome_size(),
        chriii_start = config["mat_coordinates"]["chriii_start"],
        chriii_end   = config["mat_coordinates"]["chriii_end"],
        chrv_start   = config["mat_coordinates"]["chrv_start"],
        chrv_end     = config["mat_coordinates"]["chrv_end"],
        prefix       = "{wd}/{strain}/{timepoint}_{experiment}",
    log:
        "{wd}/{strain}/logs/stage1_cov18_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        p={params.prefix}

        # BedGraph for CHRIII and CHRV loci (RPGC)
        bamCoverage -b {input.bam} -o ${{p}}_CHRIII_18nt.bedgraph \
            -of bedgraph -p {params.threads} -bs 1 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize {params.genome_size} \
            -r CHRIII:{params.chriii_start}:{params.chriii_end} >> {log} 2>&1

        bamCoverage -b {input.bam} -o ${{p}}_CHRV_18nt.bedgraph \
            -of bedgraph -p {params.threads} -bs 1 \
            --normalizeUsing RPGC \
            --effectiveGenomeSize {params.genome_size} \
            -r CHRV:{params.chrv_start}:{params.chrv_end} >> {log} 2>&1

        # Sort BedGraphs
        bedtools sort -i ${{p}}_CHRIII_18nt.bedgraph -g {input.chrom_order} > ${{p}}_CHRIII_18nt_sorted.bedgraph
        bedtools sort -i ${{p}}_CHRV_18nt.bedgraph   -g {input.chrom_order} > ${{p}}_CHRV_18nt_sorted.bedgraph

        # TSV from bedtools coverage
        bedtools coverage -a ${{p}}_CHRIII_18nt_sorted.bedgraph -b {input.bam} \
            -sorted -g {input.chrom_order} -d > ${{p}}_CHRIII_18nt.tsv 2>> {log}
        bedtools coverage -a ${{p}}_CHRV_18nt_sorted.bedgraph   -b {input.bam} \
            -sorted -g {input.chrom_order} -d > ${{p}}_CHRV_18nt.tsv 2>> {log}

        # Truncate decimals and add absolute coordinates
        awk -F'\\t' 'BEGIN{{OFS=FS}} {{split($4,a,"."); if(length(a)>1) $4=a[1]"."substr(a[2],1,5); print}}' \
            ${{p}}_CHRIII_18nt.tsv | \
            awk -F'\\t' 'BEGIN{{OFS=FS}} {{coord={params.chriii_start}+NR-1; print $1,$4,coord}}' \
            > {output.chriii_tsv}

        awk -F'\\t' 'BEGIN{{OFS=FS}} {{split($4,a,"."); if(length(a)>1) $4=a[1]"."substr(a[2],1,5); print}}' \
            ${{p}}_CHRV_18nt.tsv | \
            awk -F'\\t' 'BEGIN{{OFS=FS}} {{coord={params.chrv_start}+NR-1; print $1,$4,coord}}' \
            > {output.chrv_tsv}

        rm ${{p}}_CHRIII_18nt.bedgraph ${{p}}_CHRV_18nt.bedgraph \
           ${{p}}_CHRIII_18nt_sorted.bedgraph ${{p}}_CHRV_18nt_sorted.bedgraph \
           ${{p}}_CHRIII_18nt.tsv ${{p}}_CHRV_18nt.tsv
        """


rule coverage_18nt_nonrpgc:
    """
    Non-RPGC 18nt coverage at MAT loci for TSG/TLG/TLR timepoints.
    These TSVs are used downstream by the polymorphism diff calculation.
    Only generated for QUANT_TIMEPOINTS (TSG, TLG, TLR).
    """
    input:
        bam         = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_18nt.bam",
        bai         = "{wd}/{strain}/Alignment_data/{timepoint}_{experiment}_18nt.bai",
        chrom_order = get_chrom_order(),
    output:
        chriii_ordered = "{wd}/{strain}/MATs_quant_data/{timepoint}_{experiment}_CHRIII_18nt_nonRPGC_ordered.tsv",
        chrv_ordered   = "{wd}/{strain}/MATs_quant_data/{timepoint}_{experiment}_CHRV_18nt_nonRPGC_ordered.tsv",
    params:
        threads      = get_threads(),
        genome_size  = get_genome_size(),
        chriii_start = config["mat_coordinates"]["chriii_start"],
        chriii_end   = config["mat_coordinates"]["chriii_end"],
        chrv_start   = config["mat_coordinates"]["chrv_start"],
        chrv_end     = config["mat_coordinates"]["chrv_end"],
        prefix       = "{wd}/{strain}/{timepoint}_{experiment}",
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage1_cov18nonrpgc_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        p={params.prefix}

        bamCoverage -b {input.bam} -o ${{p}}_CHRIII_18nt_nonRPGC.bedgraph \
            -of bedgraph -p {params.threads} -bs 1 \
            --effectiveGenomeSize {params.genome_size} \
            -r CHRIII:{params.chriii_start}:{params.chriii_end} >> {log} 2>&1

        bamCoverage -b {input.bam} -o ${{p}}_CHRV_18nt_nonRPGC.bedgraph \
            -of bedgraph -p {params.threads} -bs 1 \
            --effectiveGenomeSize {params.genome_size} \
            -r CHRV:{params.chrv_start}:{params.chrv_end} >> {log} 2>&1

        bedtools sort -i ${{p}}_CHRIII_18nt_nonRPGC.bedgraph -g {input.chrom_order} \
            > ${{p}}_CHRIII_18nt_nonRPGC_sorted.bedgraph
        bedtools sort -i ${{p}}_CHRV_18nt_nonRPGC.bedgraph   -g {input.chrom_order} \
            > ${{p}}_CHRV_18nt_nonRPGC_sorted.bedgraph

        bedtools coverage -a ${{p}}_CHRIII_18nt_nonRPGC_sorted.bedgraph -b {input.bam} \
            -sorted -g {input.chrom_order} -d > ${{p}}_CHRIII_18nt_nonRPGC.tsv 2>> {log}
        bedtools coverage -a ${{p}}_CHRV_18nt_nonRPGC_sorted.bedgraph   -b {input.bam} \
            -sorted -g {input.chrom_order} -d > ${{p}}_CHRV_18nt_nonRPGC.tsv 2>> {log}

        # Truncate decimals and add absolute coordinates
        awk -F'\\t' 'BEGIN{{OFS=FS}} {{split($4,a,"."); if(length(a)>1) $4=a[1]"."substr(a[2],1,5); print}}' \
            ${{p}}_CHRIII_18nt_nonRPGC.tsv | \
            awk -F'\\t' 'BEGIN{{OFS=FS}} {{coord={params.chriii_start}+NR-1; print $1,$4,coord}}' \
            > {output.chriii_ordered}

        awk -F'\\t' 'BEGIN{{OFS=FS}} {{split($4,a,"."); if(length(a)>1) $4=a[1]"."substr(a[2],1,5); print}}' \
            ${{p}}_CHRV_18nt_nonRPGC.tsv | \
            awk -F'\\t' 'BEGIN{{OFS=FS}} {{coord={params.chrv_start}+NR-1; print $1,$4,coord}}' \
            > {output.chrv_ordered}

        rm ${{p}}_CHRIII_18nt_nonRPGC.bedgraph ${{p}}_CHRV_18nt_nonRPGC.bedgraph \
           ${{p}}_CHRIII_18nt_nonRPGC_sorted.bedgraph ${{p}}_CHRV_18nt_nonRPGC_sorted.bedgraph \
           ${{p}}_CHRIII_18nt_nonRPGC.tsv ${{p}}_CHRV_18nt_nonRPGC.tsv
        """


rule compute_polymorphism_diff:
    """
    Compute per-polymorphic-position coverage differences (mid - baseline)
    for a single timepoint × experiment.  Baseline is the mean of the
    coverage ``baseline_offset`` lines above and below the target position
    in the ordered non-RPGC TSV.

    Replaces the large awk loop in the original bash script.  Coordinates
    and offset are read from config so they are not hardcoded.
    """
    input:
        chriii = "{wd}/{strain}/MATs_quant_data/{timepoint}_{experiment}_CHRIII_18nt_nonRPGC_ordered.tsv",
        chrv   = "{wd}/{strain}/MATs_quant_data/{timepoint}_{experiment}_CHRV_18nt_nonRPGC_ordered.tsv",
    output:
        sentinel = "{wd}/{strain}/.polydiff_{timepoint}_{experiment}_done",
    params:
        chriii_coords = config["polymorphisms"]["chriii_positions"],
        chrv_coords   = config["polymorphisms"]["chrv_positions"],
        offset        = config["polymorphisms"]["baseline_offset"],
        outdir        = "{wd}/{strain}/MATs_quant_data",
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage1_polydiff_{timepoint}_{experiment}.log"
    run:
        import csv, math, logging
        from pathlib import Path

        log_path = Path(str(log[0]))
        log_path.parent.mkdir(parents=True, exist_ok=True)
        logging.basicConfig(filename=str(log_path), level=logging.INFO)

        def read_ordered_tsv(path):
            """Return list of (chr, coverage_float, coord_int) tuples."""
            rows = []
            with open(path) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    rows.append((parts[0], float(parts[1]), int(parts[2])))
            return rows

        def compute_diff(rows, target_coord, offset):
            """
            Find the row whose coordinate equals target_coord, compute
            baseline as mean of rows ±offset lines away, and return the
            diff tuple (chr, diff, coord) or None if coord not found.
            """
            coord_to_idx = {r[2]: i for i, r in enumerate(rows)}
            if target_coord not in coord_to_idx:
                return None
            mid_idx = coord_to_idx[target_coord]
            top_idx = max(0, mid_idx - offset)
            bot_idx = min(len(rows) - 1, mid_idx + offset)
            top_cov = rows[top_idx][1]
            bot_cov = rows[bot_idx][1]
            baseline = top_cov + (bot_cov - top_cov) / 2
            diff = rows[mid_idx][1] - baseline
            return (rows[mid_idx][0], diff, target_coord)

        outdir = Path(params.outdir)
        timepoint = wildcards.timepoint
        experiment = wildcards.experiment
        strain_dir = Path(wildcards.wd) / wildcards.strain

        for label, tsv_path, coords in [
            ("CHRIII", input.chriii, params.chriii_coords),
            ("CHRV",   input.chrv,   params.chrv_coords),
        ]:
            rows = read_ordered_tsv(tsv_path)
            for coord in coords:
                result = compute_diff(rows, coord, params.offset)
                if result is None:
                    logging.warning("%s coord %d not found in %s", label, coord, tsv_path)
                    continue
                chr_name, diff, c = result
                out_file = outdir / f"{timepoint}_{experiment}_{label}_diff_{coord}.tsv"
                out_file.write_text(f"{chr_name}\t{diff}\t{c}\n")
                logging.info("%s coord %d diff=%.6f", label, coord, diff)

        Path(output.sentinel).touch()


rule r_process_cov_iii_v:
    """
    Run SR_process_cov_III_V.R to compute average coverage across
    experiments at the MAT loci and produce SVG plots.
    Interface: <root_dir> <strain>
    """
    input:
        # Wait for all 75nt TSVs and ordered 18nt TSVs to be present
        tsv_75 = expand(
            "{{wd}}/{{strain}}/Coverage_data/{timepoint}_{experiment}_75nt.tsv",
            timepoint=get_timepoints(),
            experiment=get_experiments(),
        ),
        tsv_18_chriii = expand(
            "{{wd}}/{{strain}}/{timepoint}_{experiment}_CHRIII_18nt_ordered.tsv",
            timepoint=get_timepoints(),
            experiment=get_experiments(),
        ),
        tsv_18_chrv = expand(
            "{{wd}}/{{strain}}/{timepoint}_{experiment}_CHRV_18nt_ordered.tsv",
            timepoint=get_timepoints(),
            experiment=get_experiments(),
        ),
    output:
        sentinel = "{wd}/{strain}/.r_cov_iii_v_done",
    params:
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/{strain}/logs/stage1_r_cov_iii_v.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/SR_process_cov_III_V.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_process_cov_18nt:
    """
    Run SR_process_cov_III_V_18nt.R to visualise 18nt MAT quantification.
    Interface: <root_dir> <strain>
    This rule is the terminal target for the 'coverage' stage.
    """
    input:
        polydiff = expand(
            "{{wd}}/{{strain}}/.polydiff_{timepoint}_{experiment}_done",
            timepoint=QUANT_TIMEPOINTS,
            experiment=get_experiments(),
        ),
        r_cov_done = "{wd}/{strain}/.r_cov_iii_v_done",
    output:
        sentinel = "{wd}/{strain}/.stage1_done",
    params:
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/{strain}/logs/stage1_r_cov_18nt.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/SR_process_cov_III_V_18nt.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            > {log} 2>&1

        # Move plots to output subfolder
        mkdir -p {wildcards.wd}/{wildcards.strain}/Graphs_Coverage
        mv {wildcards.wd}/{wildcards.strain}/plot_*18nt*.svg \
           {wildcards.wd}/{wildcards.strain}/Graphs_Coverage/ 2>/dev/null || true

        touch {output.sentinel}
        """
