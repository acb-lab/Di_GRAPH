"""
Stage 4 — Discordant Read Mapping and Inter-chromosomal Rearrangement Analysis.

This is the most compute-intensive stage and benefits most from Snakemake
parallelism:  all strain × timepoint × experiment combinations in Group A
are independent and can run in parallel with --cores N.

Groups:
  A — 75nt paired alignment and inter-chromosomal discordant pair extraction
  B — MAT locus coverage and R visualisation (75nt)
  D — R processing: discordant pair analysis, BLAST validation, and downstream
      statistics (global distribution, categories, hotspots, radar, network)

Note: The "18nt discordant" analysis (Group C in the plan) requires a
substantial amount of additional paired-end trimming logic for 10 read-pair
subsets and is not yet implemented here; it follows the same pattern as
Group A and can be added incrementally.
"""

BLAST_OPTIONS = ["option1", "option2", "option3", "option4", "option5"]


# ---------------------------------------------------------------------------
# Group A — Alignment and SAM processing
# ---------------------------------------------------------------------------

rule trim_and_prep_disc_75nt:
    """
    Decompress paired FASTQ.gz, trim to 75nt fragments, rename read IDs
    (rA/rB suffix scheme), concatenate R1/R2 subsets, and quality-filter
    with fastp Q30.  Produces the paired filtered FASTQ files used by
    both bowtie SR and bowtie2 paired alignments.
    """
    input:
        r1 = "{wd}/{strain}/{timepoint}_{experiment}_R1.fastq.gz",
        r2 = "{wd}/{strain}/{timepoint}_{experiment}_R2.fastq.gz",
    output:
        r1_filtered = temp("{wd}/{strain}/{timepoint}_{experiment}_disc_r_R1_filtered.fastq"),
        r2_filtered = temp("{wd}/{strain}/{timepoint}_{experiment}_disc_r_R2_filtered.fastq"),
    params:
        threads  = get_threads(),
        q_thresh = config["trimming"]["quality_threshold"],
        prefix   = "{wd}/{strain}/{timepoint}_{experiment}",
    log:
        "{wd}/{strain}/logs/stage4_trim_disc_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        p={params.prefix}

        # Decompress (keep originals)
        gunzip -k {input.r1} {input.r2}
        r1_raw="$p"_R1.fastq
        r2_raw="$p"_R2.fastq

        # Trim R1 reads to two 75nt fragments
        cutadapt -j {params.threads} --cut -75 -o "$p"_R1_1.fastq "$r1_raw" >> {log} 2>&1
        cutadapt -j {params.threads} --cut  75 -o "$p"_R1_2.fastq "$r1_raw" >> {log} 2>&1
        # Trim R2 reads to two 75nt fragments
        cutadapt -j {params.threads} --cut -75 -o "$p"_R2_1.fastq "$r2_raw" >> {log} 2>&1
        cutadapt -j {params.threads} --cut  75 -o "$p"_R2_2.fastq "$r2_raw" >> {log} 2>&1
        rm "$r1_raw" "$r2_raw"

        # Rename read IDs: append rA to _1 fragments, rB to _2 fragments.
        # The sed pattern handles Novogene (LH), and other (V35/E25) sequencer headers.
        sed -E '1~4 s/^@(V35|E25|LH)[a-zA-Z0-9].*\\/1/&/; 1~4 s/\\/1/rA\\/1/' \
            "$p"_R1_1.fastq > "$p"_R1_1_A.fastq
        sed -E '1~4 s/^@(V35|E25|LH)[a-zA-Z0-9].*\\/2/&/; 1~4 s/\\/2/rA\\/2/' \
            "$p"_R2_1.fastq > "$p"_R2_1_A.fastq
        sed -E '1~4 s/^@(V35|E25|LH)[a-zA-Z0-9].*\\/1/&/; 1~4 s/\\/1/rB\\/1/' \
            "$p"_R1_2.fastq > "$p"_R1_2_B.fastq
        sed -E '1~4 s/^@(V35|E25|LH)[A-Za-z0-9].*\\/2/&/; 1~4 s/\\/2/rB\\/2/' \
            "$p"_R2_2.fastq > "$p"_R2_2_B.fastq

        rm "$p"_R1_1.fastq "$p"_R1_2.fastq "$p"_R2_1.fastq "$p"_R2_2.fastq

        # Concatenate rA and rB subsets back into paired FASTQ files
        cat "$p"_R1_1_A.fastq "$p"_R1_2_B.fastq > "$p"_disc_r_R1.fastq
        cat "$p"_R2_1_A.fastq "$p"_R2_2_B.fastq > "$p"_disc_r_R2.fastq
        rm "$p"_R1_1_A.fastq "$p"_R2_1_A.fastq "$p"_R1_2_B.fastq "$p"_R2_2_B.fastq

        # Q30 paired filter with fastp
        fastp \
            -i "$p"_disc_r_R1.fastq \
            -o {output.r1_filtered} \
            -I "$p"_disc_r_R2.fastq \
            -O {output.r2_filtered} \
            -q {params.q_thresh} -u 0 -e 30 \
            --thread {params.threads} \
            --html "$p"_disc_75nt.fastq.html \
            --json "$p"_disc_75nt.fastq.json \
            >> {log} 2>&1

        rm "$p"_disc_r_R1.fastq "$p"_disc_r_R2.fastq
        """


rule align_disc_sr_bowtie:
    """
    Independent bowtie1 SR alignment of R1 and R2 filtered reads (-m1 -v0).
    Results are used to find the subset of read names that map uniquely on
    both strands (required for concordant/discordant classification).
    """
    input:
        r1_filtered = "{wd}/{strain}/{timepoint}_{experiment}_disc_r_R1_filtered.fastq",
        r2_filtered = "{wd}/{strain}/{timepoint}_{experiment}_disc_r_R2_filtered.fastq",
        idx         = get_bowtie_index() + ".1.ebwt",
    output:
        common_reads = temp("{wd}/{strain}/{timepoint}_{experiment}_r_R12_mapped.txt"),
    params:
        threads = get_threads(),
        prefix  = get_bowtie_index(),
    log:
        "{wd}/{strain}/logs/stage4_sr_bowtie_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        p={wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}

        # Align R1 and R2 independently
        bowtie -p {params.threads} -m 1 -v 0 -S \
            {params.prefix} \
            {input.r1_filtered} > "$p"_disc_r_R1_filtered.sam 2>> {log}
        bowtie -p {params.threads} -m 1 -v 0 -S \
            {params.prefix} \
            {input.r2_filtered} > "$p"_disc_r_R2_filtered.sam 2>> {log}

        # Extract names of reads that mapped (flag -F 4)
        samtools view -F 4 "$p"_disc_r_R1_filtered.sam | cut -f1 | sed 's/\\/1//' \
            > "$p"_disc_r_R1_mapped.txt
        samtools view -F 4 "$p"_disc_r_R2_filtered.sam | cut -f1 | sed 's/\\/2//' \
            > "$p"_disc_r_R2_mapped.txt

        # Sort and find intersection (reads mapping on both R1 and R2)
        sort "$p"_disc_r_R1_mapped.txt -o "$p"_disc_r_R1_mapped_sorted.txt
        sort "$p"_disc_r_R2_mapped.txt -o "$p"_disc_r_R2_mapped_sorted.txt
        comm -12 "$p"_disc_r_R1_mapped_sorted.txt "$p"_disc_r_R2_mapped_sorted.txt \
            > {output.common_reads}

        rm "$p"_disc_r_R1_filtered.sam "$p"_disc_r_R2_filtered.sam \
           "$p"_disc_r_R1_mapped.txt "$p"_disc_r_R2_mapped.txt \
           "$p"_disc_r_R1_mapped_sorted.txt "$p"_disc_r_R2_mapped_sorted.txt
        """


rule align_disc_bowtie2_paired:
    """
    Bowtie2 paired-end alignment with strict perfect-match settings
    (--score-min C,0,0 -N 0 --end-to-end --no-1mm-upfront).
    Produces a SAM that is then processed to extract concordant and
    inter-chromosomal discordant read pairs.
    """
    input:
        r1_filtered  = "{wd}/{strain}/{timepoint}_{experiment}_disc_r_R1_filtered.fastq",
        r2_filtered  = "{wd}/{strain}/{timepoint}_{experiment}_disc_r_R2_filtered.fastq",
        common_reads = "{wd}/{strain}/{timepoint}_{experiment}_r_R12_mapped.txt",
        idx          = get_bowtie_index() + ".1.ebwt",
    output:
        concordant_count = "{wd}/{strain}/{timepoint}_{experiment}_concordant_pairs_unique_row_count.tsv",
        discordant_tsv   = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique.tsv",
    params:
        threads = get_threads(),
        prefix  = get_bowtie_index(),
    log:
        "{wd}/{strain}/logs/stage4_bt2_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        p={wildcards.wd}/{wildcards.strain}/{wildcards.timepoint}_{wildcards.experiment}

        # Bowtie2 paired alignment
        bowtie2 \
            -p {params.threads} \
            --no-1mm-upfront \
            --score-min C,0,0 \
            -N 0 \
            --end-to-end \
            --fr \
            -x {params.prefix} \
            -1 {input.r1_filtered} \
            -2 {input.r2_filtered} \
            -S "$p".sam 2>> {log}

        rm {input.r1_filtered} {input.r2_filtered}

        # Split SAM into header and alignment sections
        sed -n '1,19p' "$p".sam > "$p"_header.sam
        tail -n +20 "$p".sam  > "$p"_alignment.sam

        # ---------------------------------------------------------------
        # Concordant pair extraction:
        # Collapse quality column, pair rows, filter for proper pairs
        # ---------------------------------------------------------------
        awk -vOFS='\\t' '{{$11="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; print}}' \
            "$p"_alignment.sam | \
        cut -f1-19 | \
        sed '$!N;s/\\n/ /' | \
        awk '($3 != $7 && $7 == "=")' | \
        awk '($7 != "*")' | \
        awk '$6 == "75M"' | awk '$25 == "75M"' | \
        awk '$13 == "XN:i:0"' | awk '$32 == "XN:i:0"' | \
        sed 's/ /\\t/g' | \
        awk '($19 != "YT:Z:UP")' | awk '($38 != "YT:Z:UP")' \
        > "$p"_conc_processed.sam

        # Restore standard SAM: split into R1/R2, add header
        cut -f1-19 "$p"_conc_processed.sam > "$p"_conc_R1.sam
        cut -f20-38 "$p"_conc_processed.sam > "$p"_conc_R2.sam
        cat "$p"_conc_R1.sam "$p"_conc_R2.sam > "$p"_conc_R12.sam
        cat "$p"_header.sam "$p"_conc_R12.sam > "$p"_concordant_pairs.sam
        rm "$p"_conc_processed.sam "$p"_conc_R1.sam "$p"_conc_R2.sam "$p"_conc_R12.sam

        # Filter to reads that also mapped uniquely in SR mode
        samtools view -N {input.common_reads} \
            -o "$p"_conc_R12_R1.sam "$p"_concordant_pairs.sam
        samtools view -N {input.common_reads} \
            -o "$p"_conc_R12_R2.sam "$p"_concordant_pairs.sam
        cat "$p"_conc_R12_R1.sam "$p"_conc_R12_R2.sam > "$p"_conc_R12_final.sam
        cat "$p"_header.sam "$p"_conc_R12_final.sam > "$p"_concordant_pairs_full.sam

        # Deduplicate
        awk '!seen[$0]++' "$p"_concordant_pairs_full.sam > "$p"_concordant_pairs_unique.sam
        # Count (excluding header)
        tail -n +20 "$p"_concordant_pairs_unique.sam | wc -l | awk '{{print $1}}' \
            > {output.concordant_count}

        # Clean up concordant intermediates
        rm "$p"_concordant_pairs.sam "$p"_conc_R12_R1.sam "$p"_conc_R12_R2.sam \
           "$p"_conc_R12_final.sam "$p"_concordant_pairs_full.sam \
           "$p"_concordant_pairs_unique.sam

        # ---------------------------------------------------------------
        # Inter-chromosomal discordant pair extraction:
        # $3 != $7 (different chromosomes), $7 != "=" and $7 != "*"
        # ---------------------------------------------------------------
        awk -vOFS='\\t' '{{$11="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; print}}' \
            "$p"_alignment.sam | \
        cut -f1-19 | \
        sed '$!N;s/\\n/ /' | \
        awk '($3 != $7 && $7 != "=")' | \
        awk '($7 != "*")' | \
        awk '$6 == "75M"' | awk '$25 == "75M"' | \
        awk '$13 == "XN:i:0"' | awk '$32 == "XN:i:0"' | \
        sed 's/ /\\t/g' | \
        awk '($19 != "YT:Z:UP")' | awk '($38 != "YT:Z:UP")' \
        > "$p"_disc_processed.sam

        cut -f1-19 "$p"_disc_processed.sam > "$p"_disc_R1.sam
        cut -f20-38 "$p"_disc_processed.sam > "$p"_disc_R2.sam
        cat "$p"_disc_R1.sam "$p"_disc_R2.sam > "$p"_disc_R12.sam
        cat "$p"_header.sam "$p"_disc_R12.sam > "$p"_inter_disc_pairs.sam
        rm "$p"_disc_processed.sam "$p"_disc_R1.sam "$p"_disc_R2.sam "$p"_disc_R12.sam

        samtools view -N {input.common_reads} \
            -o "$p"_disc_R12_R1.sam "$p"_inter_disc_pairs.sam
        samtools view -N {input.common_reads} \
            -o "$p"_disc_R12_R2.sam "$p"_inter_disc_pairs.sam
        cat "$p"_disc_R12_R1.sam "$p"_disc_R12_R2.sam > "$p"_disc_R12_final.sam
        cat "$p"_header.sam "$p"_disc_R12_final.sam > "$p"_inter_disc_full.sam

        awk '!seen[$0]++' "$p"_inter_disc_full.sam > "$p"_inter_discordant_pairs_unique.sam

        # Extract columns: read_name, chr_A, pos_A, chr_B, pos_B, sequence
        tail -n +20 "$p"_inter_discordant_pairs_unique.sam | \
        awk 'BEGIN{{OFS="\\t"}} {{print $1,$3,$4,$7,$8,$10}}' \
        > {output.discordant_tsv}

        # Clean up
        rm "$p"_alignment.sam "$p"_header.sam "$p"_inter_disc_pairs.sam \
           "$p"_disc_R12_R1.sam "$p"_disc_R12_R2.sam "$p"_disc_R12_final.sam \
           "$p"_inter_disc_full.sam "$p"_inter_discordant_pairs_unique.sam \
           "$p"_inter_disc_pairs.sam "$p".sam {input.common_reads} 2>/dev/null || true
        """


# ---------------------------------------------------------------------------
# Group B — MAT locus analysis (75nt discordant)
# ---------------------------------------------------------------------------

rule r_plot_mat_coverage_75nt:
    """
    Plot MAT locus coverage from 75nt discordant reads.
    Interface: <root_dir> <strain>
    """
    input:
        tsv_files = expand(
            "{{wd}}/{{strain}}/{timepoint}_{experiment}_inter_discordant_pairs_unique.tsv",
            timepoint=get_timepoints(),
            experiment=get_experiments(),
        ),
    output:
        sentinel = "{wd}/{strain}/.disc_mat_cov_75nt_done",
    params:
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/{strain}/logs/stage4_mat_cov_75nt.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_plot_MATa_cov.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_plot_mat_coverage_18nt:
    """
    Plot MAT locus coverage from 18nt discordant reads.
    Interface: <root_dir> <strain>
    """
    input:
        mat_done = "{wd}/{strain}/.disc_mat_cov_75nt_done",
    output:
        sentinel = "{wd}/{strain}/.disc_mat_cov_18nt_done",
    params:
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/{strain}/logs/stage4_mat_cov_18nt.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_plot_MATa_cov_r18.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            > {log} 2>&1
        touch {output.sentinel}
        """


# ---------------------------------------------------------------------------
# Group D — R processing pipeline (discordant pair analysis and statistics)
# ---------------------------------------------------------------------------

rule r_process_inter_discordant_10kb:
    """
    Filter discordant pairs for ≥10kb inter-chromosomal distance and
    associate to genomic features.
    Interface: <strain> <sample> <experiment> <file_T0> <file_sample> <category_path>
    """
    input:
        t0_tsv  = "{wd}/{strain}/T0_{experiment}_inter_discordant_pairs_unique.tsv",
        smp_tsv = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique.tsv",
    output:
        processed_tsv = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed.tsv",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_process_disc_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_process_inter_discordant_pairs_10kb.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.timepoint} \
            {wildcards.experiment} \
            {input.t0_tsv} \
            {input.smp_tsv} \
            {params.category_path} \
            > {log} 2>&1
        """


rule blast_cross_validate:
    """
    BLAST cross-validation of inter-chromosomal discordant pairs.
    For each read pair, query its sequence against the corresponding
    chromosome FASTA to discard PCR artefacts and sequencing errors.
    The loop logic is kept verbatim from the original bash script.
    """
    input:
        blast_tsv = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed_blast_chromosome.tsv",
    output:
        combined = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed_blast_chromosome_combined_results.tsv",
    params:
        blast_dir = get_blast_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_blast_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        cp {input.blast_tsv} {params.blast_dir}/
        file_name=$(basename "{input.blast_tsv}" .tsv)
        output_dir={params.blast_dir}/1.blast_results
        mkdir -p "$output_dir"

        (
            cd {params.blast_dir} || exit

            while IFS=$'\\t' read -r nombre_read query_A reference_B || [ -n "$nombre_read" ]; do
                [ -z "$nombre_read" ] || [ -z "$query_A" ] || [ -z "$reference_B" ] && continue
                reference_ext="_v9.fasta"
                reference_file_path="${{reference_B}}${{reference_ext}}"
                [ -f "$reference_file_path" ] || continue

                random_id=$RANDOM
                echo -e ">${{nombre_read}}\\n${{query_A}}" > "${{file_name}}".fasta
                makeblastdb -in "$reference_file_path" -dbtype nucl -out "${{reference_B}}" \
                    >> {log} 2>&1
                blastn -db "${{reference_B}}" -query "${{file_name}}".fasta \
                    -dust no \
                    -outfmt "7 qseqid sseqid mismatch length sstart send qseq sseq" \
                    -out "${{nombre_read}}_${{reference_B}}_randomid${{random_id}}_results.tsv" \
                    >> {log} 2>&1
                mv *_results.tsv "$output_dir"
                rm "${{file_name}}".fasta
                rm -f *.nhr *.nin *.nsq *.ndb *.njs *.not *.ntf *.nto

                result_file="${{output_dir}}/${{nombre_read}}_${{reference_B}}_randomid${{random_id}}_results.tsv"
                file_base="${{result_file%.tsv}}"
                grep -E 'Query|Database|hits|E25|V35' "$result_file" \
                    > "${{file_base}}_selected.tsv"
                tr ' ' '\\t' < "${{file_base}}_selected.tsv" \
                    > "${{file_base}}_selected_tab.tsv"
                tr '\\n ' '\\t' < "${{file_base}}_selected_tab.tsv" \
                    > "${{file_base}}_selected_tab_n.tsv"
                cut -f3,6,8,13- "${{file_base}}_selected_tab_n.tsv" \
                    > "${{file_base}}_results_processed.tsv"
                rm -f "${{file_base}}_selected.tsv" "${{file_base}}_selected_tab.tsv" \
                      "${{file_base}}_selected_tab_n.tsv" "$result_file"
            done < {input.blast_tsv}

            # Combine all processed results
            find "$output_dir" -name '*_results_results_processed.tsv' -print0 | \
                xargs -0 cat > {output.combined}
            find "$output_dir" -name '*_results_results_processed.tsv' -delete
            rm -r "$output_dir"
            rm -f "$(basename {input.blast_tsv})"
        )
        """


rule r_validate_inter_discordant_10kb:
    """
    Validate discordant pairs against BLAST results and category annotations.
    Interface: <strain> <sample> <experiment> <blast_results> <orig_tsv> <control_tsv> <category_path>
    """
    input:
        blast_combined = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed_blast_chromosome_combined_results.tsv",
        orig_tsv       = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed.tsv",
        control_tsv    = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed_control.tsv",
    output:
        sentinel = "{wd}/{strain}/.disc_validated_{timepoint}_{experiment}_done",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_validate_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_validate_inter_discordant_pairs_10kb.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.timepoint} \
            {wildcards.experiment} \
            {input.blast_combined} \
            {input.orig_tsv} \
            {input.control_tsv} \
            {params.category_path} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_validate_blast_levels:
    """
    Classify validated discordant pairs into 5 stringency levels.
    Interface: <strain> <sample> <experiment> <category_path> <blast_results> <orig_tsv> <control_tsv>
    """
    input:
        validated_done = "{wd}/{strain}/.disc_validated_{timepoint}_{experiment}_done",
        blast_combined = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed_blast_chromosome_combined_results.tsv",
        orig_tsv       = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed.tsv",
        control_tsv    = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed_control.tsv",
    output:
        sentinel = "{wd}/{strain}/.disc_levels_{timepoint}_{experiment}_done",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_levels_{timepoint}_{experiment}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_validate_inter_discordant_pairs_levels_1feature_10kb.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.timepoint} \
            {wildcards.experiment} \
            {params.category_path} \
            {input.blast_combined} \
            {input.orig_tsv} \
            {input.control_tsv} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_process_valid_levels:
    """
    Per blast_option processing of validated discordant pairs.
    Interface: <strain> <sample> <experiment> <blast_option> <category_path> <file1> <file2> <control_tsv>
    """
    input:
        levels_done = "{wd}/{strain}/.disc_levels_{timepoint}_{experiment}_done",
        file1       = "{wd}/{strain}/{timepoint}_{experiment}_{blast_option}_inter_discordant_pairs_unique_processed_valid.tsv",
        file2       = "{wd}/{strain}/{timepoint}_{experiment}_{blast_option}_inter_discordant_pairs_unique_processed_valid_error_rate.tsv",
        control_tsv = "{wd}/{strain}/{timepoint}_{experiment}_inter_discordant_pairs_unique_processed_control.tsv",
    output:
        sentinel = "{wd}/{strain}/.disc_proc_valid_{timepoint}_{experiment}_{blast_option}_done",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    wildcard_constraints:
        timepoint    = "TSG|TLG|TLR",
        blast_option = "option[1-5]",
    log:
        "{wd}/{strain}/logs/stage4_proc_valid_{timepoint}_{experiment}_{blast_option}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_process_valid_inter_discordant_pairs_levels.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.timepoint} \
            {wildcards.experiment} \
            {wildcards.blast_option} \
            {params.category_path} \
            {input.file1} \
            {input.file2} \
            {input.control_tsv} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_global_distribution:
    """Calculate global distribution of discordant pairs. Interface: <root_dir> <strain> <sample>"""
    input:
        proc_done = expand(
            "{{wd}}/{{strain}}/.disc_validated_{{timepoint}}_{experiment}_done",
            experiment=get_experiments(),
        ),
    output:
        sentinel = "{wd}/{strain}/.disc_global_dist_{timepoint}_done",
    params:
        scripts_dir = get_scripts_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_global_dist_{timepoint}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_global_distribution.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.timepoint} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_global_distribution_levels:
    """Global distribution per blast_option. Interface: <root_dir> <strain> <sample> <blast_option>"""
    input:
        all_valid_done = expand(
            "{{wd}}/{{strain}}/.disc_proc_valid_{{timepoint}}_{experiment}_{{blast_option}}_done",
            experiment=get_experiments(),
        ),
    output:
        sentinel = "{wd}/{strain}/.disc_global_dist_levels_{timepoint}_{blast_option}_done",
    params:
        scripts_dir = get_scripts_dir(),
    wildcard_constraints:
        timepoint    = "TSG|TLG|TLR",
        blast_option = "option[1-5]",
    log:
        "{wd}/{strain}/logs/stage4_global_dist_levels_{timepoint}_{blast_option}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_global_distribution_levels.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.timepoint} \
            {wildcards.blast_option} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_category_distribution:
    """Category distribution. Interface: <root_dir> <strain> <sample>"""
    input:
        global_done = "{wd}/{strain}/.disc_global_dist_{timepoint}_done",
    output:
        sentinel = "{wd}/{strain}/.disc_cat_dist_{timepoint}_done",
    params:
        scripts_dir = get_scripts_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_cat_dist_{timepoint}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_category_distribution.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.timepoint} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_matrix_all_exp:
    """Merge all experiments into a rearrangement matrix. Interface: <root_dir> <strain> <sample> <category_path>"""
    input:
        cat_done = "{wd}/{strain}/.disc_cat_dist_{timepoint}_done",
    output:
        sentinel = "{wd}/{strain}/.disc_matrix_{timepoint}_done",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_matrix_{timepoint}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_matrix_all_exp.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.timepoint} \
            {params.category_path} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_matrix_all_exp_levels:
    """Matrix per blast_option. Interface: <root_dir> <strain> <sample> <blast_option> <category_path>"""
    input:
        levels_done = "{wd}/{strain}/.disc_global_dist_levels_{timepoint}_{blast_option}_done",
    output:
        sentinel = "{wd}/{strain}/.disc_matrix_levels_{timepoint}_{blast_option}_done",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    wildcard_constraints:
        timepoint    = "TSG|TLG|TLR",
        blast_option = "option[1-5]",
    log:
        "{wd}/{strain}/logs/stage4_matrix_levels_{timepoint}_{blast_option}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_matrix_all_exp_levels.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.timepoint} \
            {wildcards.blast_option} \
            {params.category_path} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_distrib_category_a:
    """Category_A distribution. Interface: <root_dir> <strain>"""
    input:
        matrices_done = expand(
            "{{wd}}/{{strain}}/.disc_matrix_{timepoint}_done",
            timepoint=QUANT_TIMEPOINTS,
        ),
    output:
        sentinel = "{wd}/{strain}/.disc_distrib_catA_done",
    params:
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/{strain}/logs/stage4_distrib_catA.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_distrib_categoryA.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_identify_hotspots:
    """Identify recombination hotspots. Interface: <root_dir> <strain> <sample> <path_to_freqs>"""
    input:
        matrix_done = "{wd}/{strain}/.disc_matrix_{timepoint}_done",
        freqs       = "{wd}/{strain}/{timepoint}_inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments.tsv",
    output:
        sentinel = "{wd}/{strain}/.disc_hotspots_{timepoint}_done",
    params:
        scripts_dir = get_scripts_dir(),
    wildcard_constraints:
        timepoint = "TSG|TLG|TLR",
    log:
        "{wd}/{strain}/logs/stage4_hotspots_{timepoint}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_identify_hotspots.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.timepoint} \
            {input.freqs} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_identify_hotspots_levels:
    """Hotspots per blast_option. Interface: <root_dir> <strain> <sample> <blast_option> <path_to_freqs>"""
    input:
        matrix_done = "{wd}/{strain}/.disc_matrix_levels_{timepoint}_{blast_option}_done",
        freqs       = "{wd}/{strain}/{timepoint}_{blast_option}_inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments.tsv",
    output:
        sentinel = "{wd}/{strain}/.disc_hotspots_levels_{timepoint}_{blast_option}_done",
    params:
        scripts_dir = get_scripts_dir(),
    wildcard_constraints:
        timepoint    = "TSG|TLG|TLR",
        blast_option = "option[1-5]",
    log:
        "{wd}/{strain}/logs/stage4_hotspots_levels_{timepoint}_{blast_option}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_identify_hotspots_levels.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.timepoint} \
            {wildcards.blast_option} \
            {input.freqs} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_hotspots_distribution_levels:
    """
    Hotspot distribution across timepoints (whole WD).
    Interface: <root_dir>
    """
    input:
        all_hotspot_levels = expand(
            "{wd}/{strain}/.disc_hotspots_levels_{timepoint}_{blast_option}_done",
            wd=get_working_dir(),
            strain=get_sample_names(),
            timepoint=QUANT_TIMEPOINTS,
            blast_option=BLAST_OPTIONS,
        ),
    output:
        sentinel = f"{get_working_dir()}/.disc_hotspots_dist_done",
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    log:
        f"{get_working_dir()}/logs/stage4_hotspots_dist.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_hotspots_distribution_levels.R \
            {params.wd} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_recombination_rate:
    """
    Calculate genome-wide recombination rate across all strains.
    Interface: <root_dir> <reference_strain>
    """
    input:
        hotspots_dist_done = f"{get_working_dir()}/.disc_hotspots_dist_done",
    output:
        sentinel = f"{get_working_dir()}/.disc_recomb_rate_done",
    params:
        scripts_dir      = get_scripts_dir(),
        wd               = get_working_dir(),
        reference_strain = lambda wc: (
            next((s["name"] for s in config["samples"] if s.get("is_reference")), "no")
        ),
    log:
        f"{get_working_dir()}/logs/stage4_recomb_rate.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_recombination_rate.R \
            {params.wd} \
            {params.reference_strain} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_recombination_rate_conc_disc:
    """
    Per-strain concordant/discordant recombination rate.
    Interface: <root_dir> <strain> <wd_dir>
    """
    input:
        rate_done = f"{get_working_dir()}/.disc_recomb_rate_done",
    output:
        sentinel = "{wd}/{strain}/.disc_recomb_rate_conc_disc_done",
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    log:
        "{wd}/{strain}/logs/stage4_recomb_rate_conc_disc.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_recombination_rate_conc_disc.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {params.wd} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_distrib_read_number_category_a:
    """
    Category_A read number distribution.
    Interface: <root_dir> <strain> <wd_dir>
    """
    input:
        distrib_done = "{wd}/{strain}/.disc_distrib_catA_done",
    output:
        sentinel = "{wd}/{strain}/.disc_read_num_catA_done",
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    log:
        "{wd}/{strain}/logs/stage4_read_num_catA.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_distrib_read_number_categoryA.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {params.wd} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_error_rate:
    """
    Estimate false-positive (error) rate.
    Interface: <root_dir> <strain> <wd_dir>
    """
    input:
        rate_conc_disc_done = "{wd}/{strain}/.disc_recomb_rate_conc_disc_done",
    output:
        sentinel = "{wd}/{strain}/.disc_error_rate_done",
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    log:
        "{wd}/{strain}/logs/stage4_error_rate.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_error_rate.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {params.wd} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_radar_perc_levels:
    """
    Radar plot of category percentages per blast_option.
    Interface: <root_dir> <strain> <wd_dir> <blast_option>
    """
    input:
        error_done = "{wd}/{strain}/.disc_error_rate_done",
    output:
        sentinel = "{wd}/{strain}/.disc_radar_perc_{blast_option}_done",
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    wildcard_constraints:
        blast_option = "option[1-5]",
    log:
        "{wd}/{strain}/logs/stage4_radar_perc_{blast_option}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_radar_perc_levels.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {params.wd} \
            {wildcards.blast_option} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_radar_number_levels:
    """
    Radar plot of read counts.
    Interface: <root_dir> <strain> <wd_dir>
    """
    input:
        read_num_done = "{wd}/{strain}/.disc_read_num_catA_done",
    output:
        sentinel = "{wd}/{strain}/.disc_radar_number_done",
    params:
        scripts_dir = get_scripts_dir(),
        wd          = get_working_dir(),
    log:
        "{wd}/{strain}/logs/stage4_radar_number.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/DISC_calculate_radar_number_levels.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {params.wd} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_discordant_network:
    """
    Build discordant read interaction network.
    Interface: <root_dir> <strain> <sample>
    Terminal rule for stage 4 per strain.
    """
    input:
        radar_perc_done = expand(
            "{{wd}}/{{strain}}/.disc_radar_perc_{blast_option}_done",
            blast_option=BLAST_OPTIONS,
        ),
        radar_number_done = "{wd}/{strain}/.disc_radar_number_done",
        hotspots_done = expand(
            "{{wd}}/{{strain}}/.disc_hotspots_{timepoint}_done",
            timepoint=QUANT_TIMEPOINTS,
        ),
        mat_cov_done = "{wd}/{strain}/.disc_mat_cov_18nt_done",
    output:
        sentinel = "{wd}/{strain}/.stage4_done",
    params:
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/{strain}/logs/stage4_network.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        for sample in TSG TLG TLR; do
            Rscript {params.scripts_dir}/R/DISC_calculate_discordant_network.R \
                {wildcards.wd}/{wildcards.strain}/ \
                {wildcards.strain} \
                "$sample" >> {log} 2>&1
        done

        # Organise outputs
        mkdir -p {wildcards.wd}/{wildcards.strain}/Discordant_global_analysis/Global_Data_75nt
        mkdir -p {wildcards.wd}/{wildcards.strain}/Discordant_global_analysis/Global_Plots_75nt
        mv {wildcards.wd}/{wildcards.strain}/*.svg \
           {wildcards.wd}/{wildcards.strain}/Discordant_global_analysis/Global_Plots_75nt/ 2>/dev/null || true

        touch {output.sentinel}
        """
