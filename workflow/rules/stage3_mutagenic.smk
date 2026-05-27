"""
Stage 3 — Mutagenic Rate at the HO Site.

For each strain × timepoint (T0/TLG/TLR) × experiment:
  1. Build BWA index on CHRIII FASTA (once).
  2. Decompress → concatenate R1+R2 → fastp Q30 filter.
  3. BWA mem alignment → filter CHRIII reads → sorted BAM.
  4. Extract reads spanning the HO site (CHRIII:200753).
  5. Python scripts classify reads by nucleotide at 200689 and 200753.
  6. igvtools + bcftools assess deletion/insertion and nucleotide variants.
  7. R scripts compute repair pathway frequencies and generate plots.

R script interfaces:
  MUT_calc_rep_path_freq.R      <root_dir> <strain> <wd_dir>
  MUT_plot_rep_path_comparison.R <root_dir>

Python script interfaces:
  MUT_extract_200689.py  <input.bam> <output_A.bam> <output_G.bam>
  MUT_extract_200753.py  <input.bam> <output_T.bam> <output_C.bam>
"""


rule build_bwa_index:
    """Build the BWA index on the CHRIII reference FASTA (done once)."""
    input:
        fasta = get_genome_fasta_chriii(),
    output:
        # BWA produces 5 index files; track .bwt as proxy
        idx = get_genome_fasta_chriii() + ".bwt",
    log:
        f"{get_working_dir()}/logs/build_bwa_index.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        "bwa index {input.fasta} > {log} 2>&1"


rule filter_and_align_bwa:
    """
    Decompress + concatenate R1/R2 → fastp Q30 filter → BWA mem alignment
    → extract CHRIII-only reads → sorted BAM.
    Runs for T0, TLG, TLR only.
    """
    input:
        r1  = "{wd}/{strain}/{timepoint}_{experiment}_R1.fastq.gz",
        r2  = "{wd}/{strain}/{timepoint}_{experiment}_R2.fastq.gz",
        idx = get_genome_fasta_chriii() + ".bwt",
    output:
        bam = "{wd}/{strain}/{timepoint}_{experiment}_processed.bam",
        bai = "{wd}/{strain}/{timepoint}_{experiment}_processed.bai",
    params:
        threads   = get_threads(),
        ref_chriii = get_genome_fasta_chriii(),
        q_thresh  = config["trimming"]["quality_threshold"],
    wildcard_constraints:
        timepoint = "T0|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage3_bwa_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        export LC_ALL=C LANG=C

        # Decompress, concatenate
        gunzip -k {input.r1} {input.r2}
        r1={wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_R1.fastq
        r2={wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_R2.fastq
        cat "$r1" "$r2" > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}.fastq
        rm "$r1" "$r2"

        # Q30 filter
        fastp \
            -i {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}.fastq \
            -o {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.fastq \
            -q {params.q_thresh} -u 0 -e 30 \
            --thread {params.threads} \
            --html {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.fastq.html \
            --json {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.fastq.json \
            >> {log} 2>&1
        rm {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}.fastq

        # BWA alignment
        bwa mem -t {params.threads} {params.ref_chriii} \
            {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.fastq \
            > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.sam 2>> {log}

        # Extract header (first 2 lines) and CHRIII reads
        sam={wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.sam
        sed -n '1,2p' "$sam" > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_cabecera.sam
        awk '$3=="CHRIII"' "$sam" > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30_CHRIII.sam
        cat {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_cabecera.sam \
            {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30_CHRIII.sam \
            > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30_CHRIII_header.sam

        # Sort and index BAM
        samtools sort -@ {params.threads} \
            -o {output.bam} \
            {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30_CHRIII_header.sam 2>> {log}
        samtools index {output.bam} {output.bai}

        # Clean up intermediates
        rm {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.fastq \
           {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30.sam \
           {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30_CHRIII.sam \
           {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_cabecera.sam \
           {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_Q30_CHRIII_header.sam
        """


rule extract_ho_site_reads:
    """Extract reads overlapping the HO cut site (CHRIII:200753)."""
    input:
        bam = "{wd}/{strain}/{timepoint}_{experiment}_processed.bam",
        bai = "{wd}/{strain}/{timepoint}_{experiment}_processed.bai",
    output:
        bam = "{wd}/{strain}/{timepoint}_{experiment}_HOs.bam",
        bai = "{wd}/{strain}/{timepoint}_{experiment}_HOs.bai",
    params:
        threads  = get_threads(),
        ho_site  = config["mat_coordinates"]["ho_site"],
    wildcard_constraints:
        timepoint = "T0|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage3_ho_extract_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        samtools view \
            {input.bam} "CHRIII:{params.ho_site}-{params.ho_site}" \
            -h -O SAM \
            > {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_HOs.sam 2>> {log}

        samtools sort -@ {params.threads} \
            -o {output.bam} \
            {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_HOs.sam 2>> {log}
        samtools index {output.bam} {output.bai}

        rm {wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}_HOs.sam
        """


rule python_extract_200689:
    """
    Classify HO-site reads by nucleotide at CHRIII position 200689 (A or G).
    Uses MUT_extract_200689.py via pysam.
    """
    input:
        bam = "{wd}/{strain}/{timepoint}_{experiment}_HOs.bam",
    output:
        bam_A = "{wd}/{strain}/{timepoint}_{experiment}_200689_A.bam",
        bam_G = "{wd}/{strain}/{timepoint}_{experiment}_200689_G.bam",
    params:
        script = f"{get_scripts_dir()}/python/MUT_extract_200689.py",
    wildcard_constraints:
        timepoint = "T0|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage3_extract200689_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        python3 {params.script} \
            {input.bam} {output.bam_A} {output.bam_G} \
            > {log} 2>&1
        samtools index {output.bam_A}
        samtools index {output.bam_G}
        """


rule python_extract_200753:
    """
    From reads with A at 200689, classify by nucleotide at 200753 (T or C).
    Uses MUT_extract_200753.py via pysam.
    """
    input:
        bam = "{wd}/{strain}/{timepoint}_{experiment}_200689_A.bam",
    output:
        bam_T = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_200753_T.bam",
        bam_C = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_200753_C.bam",
    params:
        script = f"{get_scripts_dir()}/python/MUT_extract_200753.py",
    wildcard_constraints:
        timepoint = "T0|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage3_extract200753_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        python3 {params.script} \
            {input.bam} {output.bam_T} {output.bam_C} \
            > {log} 2>&1
        samtools index {output.bam_T}
        samtools index {output.bam_C}
        """


rule merge_bam_variants:
    """
    Merge classified BAM files to produce combined variant groups
    used by igvtools and bcftools.
    """
    input:
        bam_A   = "{wd}/{strain}/{timepoint}_{experiment}_200689_A.bam",
        bam_T   = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_200753_T.bam",
        bam_C   = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_200753_C.bam",
        bam_G   = "{wd}/{strain}/{timepoint}_{experiment}_200689_G.bam",
    output:
        bam_A_CT   = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_CT.bam",
        bam_A_CT_G = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_CT_G.bam",
    wildcard_constraints:
        timepoint = "T0|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage3_merge_bam_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        # A_CT = reads with A at 200689 (T+C at 200753 combined)
        samtools merge -f {output.bam_A_CT} {input.bam_T} {input.bam_C} >> {log} 2>&1
        samtools index {output.bam_A_CT}

        # A_CT_G = all reads (A_CT + G at 200689)
        samtools merge -f {output.bam_A_CT_G} {output.bam_A_CT} {input.bam_G} >> {log} 2>&1
        samtools index {output.bam_A_CT_G}
        """


rule igvtools_and_bcftools:
    """
    For each variant class, run igvtools count (WIG) and bcftools
    mpileup/call to detect nucleotide variants at the HO site region.
    Produces a BCF TSV per variant class.
    """
    input:
        bam_A_CT   = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_CT.bam",
        bam_A_CT_G = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_CT_G.bam",
        bam_A_T    = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_200753_T.bam",
        bam_A_C    = "{wd}/{strain}/{timepoint}_{experiment}_200689_A_200753_C.bam",
        bam_G      = "{wd}/{strain}/{timepoint}_{experiment}_200689_G.bam",
    output:
        sentinel = "{wd}/{strain}/.igv_bcf_{timepoint}_{experiment}_done",
    params:
        ref_chriii   = get_genome_fasta_chriii(),
        ho_region    = "CHRIII:{chriii_start}-{chriii_end}".format(
            chriii_start=config["mat_coordinates"]["ho_site"] - 10,
            chriii_end  =config["mat_coordinates"]["ho_site"] + 10,
        ),
        prefix       = "{wd}/{strain}/{timepoint}_{experiment}",
    wildcard_constraints:
        timepoint = "T0|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage3_igv_bcf_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        p={params.prefix}
        ref={params.ref_chriii}

        # Map variant class → BAM file
        declare -A bam_map
        bam_map["A_CT"]={input.bam_A_CT}
        bam_map["A_CT_G"]={input.bam_A_CT_G}
        bam_map["A_200753_T"]={input.bam_A_T}
        bam_map["A_200753_C"]={input.bam_A_C}
        bam_map["G"]={input.bam_G}

        for poly in A_CT A_CT_G A_200753_T A_200753_C G; do
            bam="${{bam_map[$poly]}}"
            wig="$p"_200689_${{poly}}.wig
            tsv="$p"_${{poly}}_BCF.tsv

            # igvtools WIG (1bp resolution over HO site)
            igvtools count -w 1 --bases \
                --query CHRIII:{params.ho_region} \
                "$bam" "$wig" "$ref" >> {log} 2>&1

            # BCFtools variant calling
            vcf="$p"_${{poly}}_variants.vcf.gz
            bcftools mpileup -Ou -f "$ref" \
                -r CHRIII:{params.ho_region} "$bam" | \
            bcftools call -mA --ploidy 1 -Oz -o "$vcf" >> {log} 2>&1
            bcftools index "$vcf"

            vcf_af="$p"_${{poly}}_with_AF.vcf.gz
            bcftools +fill-tags "$vcf" -Oz -o "$vcf_af" -- -t AF,AN,AC >> {log} 2>&1
            bcftools index "$vcf_af"

            # Extract numeric summary to TSV
            bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/DP\\t%INFO/DP4\\n' "$vcf" | \
            awk 'BEGIN{{OFS="\\t"}} {{
                split($6,dp4,",");
                ref_count=dp4[1]+dp4[2];
                alt_count=dp4[3]+dp4[4];
                total=ref_count+alt_count;
                af=(total==0)?0:int((alt_count/total)*100);
                print $1,$2,$3,$4,$5,af
            }}' > "$tsv"

            rm -f "$vcf" "$vcf_af" "${{vcf%.vcf.gz}}.csi" "${{vcf_af%.vcf.gz}}.csi"
            rm -f "$wig"
        done

        touch {output.sentinel}
        """


rule r_calc_repair_pathway:
    """
    Calculate repair pathway frequencies per strain.
    Interface: <root_dir> <strain> <wd_dir>
    """
    input:
        igv_done = expand(
            "{{wd}}/{{strain}}/.igv_bcf_{timepoint}_{experiment}_done",
            timepoint=MUT_TIMEPOINTS,
            experiment=get_experiments(),
        ),
    output:
        sentinel = "{wd}/{strain}/.stage3_repair_freq_done",
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    log:
        "{wd}/{strain}/logs/stage3_repair_freq.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/MUT_calc_rep_path_freq.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {params.wd} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_plot_repair_comparison:
    """
    Plot repair pathway comparison across all strains in the working directory.
    Interface: <root_dir>
    Terminal rule for stage 3 — produces .stage3_done sentinel per strain.
    """
    input:
        all_freq_done = expand(
            "{wd}/{strain}/.stage3_repair_freq_done",
            wd=get_working_dir(),
            strain=get_sample_names(),
        ),
    output:
        expand(
            "{wd}/{strain}/.stage3_done",
            wd=get_working_dir(),
            strain=get_sample_names(),
        ),
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    log:
        f"{get_working_dir()}/logs/stage3_repair_comparison.log"
    conda:
        "../envs/digraph.yaml"
    run:
        import subprocess, pathlib
        r = subprocess.run(
            [
                "Rscript",
                f"{params.scripts_dir}/R/MUT_plot_rep_path_comparison.R",
                params.wd,
            ],
            capture_output=True, text=True,
        )
        pathlib.Path(log[0]).parent.mkdir(parents=True, exist_ok=True)
        pathlib.Path(log[0]).write_text(r.stdout + r.stderr)
        if r.returncode != 0:
            raise RuntimeError(f"MUT_plot_rep_path_comparison.R failed: {r.stderr[:500]}")

        # Organise outputs and write sentinels
        wd = pathlib.Path(params.wd)
        for strain in get_sample_names():
            strain_dir = wd / strain
            ho_dir = strain_dir / "HO_mut_rate"
            for sub in ("FASTQP", "SVG", "TSV"):
                (ho_dir / sub).mkdir(parents=True, exist_ok=True)
            for ext, sub in [("*.fastq.html", "FASTQP"), ("*.fastq.json", "FASTQP"),
                              ("*.svg", "SVG"), ("*.tsv", "TSV")]:
                for f in strain_dir.glob(ext):
                    f.rename(ho_dir / sub / f.name)
            (strain_dir / ".stage3_done").touch()
