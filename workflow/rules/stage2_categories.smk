"""
Stage 2 — Genomic Categories Analysis.

For each strain, runs four R scripts that profile coverage changes across
the 13 genomic feature categories (ORF, LTR, TEG, Ty, tRNA, rRNA, ncRNA,
snRNA, snoRNA, ARS, Cen, Tel, Int) at each pairwise timepoint comparison
(T0vsTSG, T0vsTLG, T0vsTLR).

R script interfaces (verified from source):
  SR_process_categories.R  <root_dir> <strain> <category> <category_path>
  SR_order_categories.R    <root_dir> <strain> <category_path> <analysis_suffix>
  SR_plot_categories.R     <root_dir> <strain> <analysis_suffix>
  SR_plot_gal_vs_raf.R     <root_dir> <strain>
"""

ANALYSIS_SUFFIXES = ["T0vsTSG", "T0vsTLG", "T0vsTLR"]


rule r_process_categories:
    """
    Compute per-category average coverage fingerprints for a single category.
    One rule instance per (strain, category).
    """
    input:
        stage1_done = "{wd}/{strain}/.stage1_done",
    output:
        sentinel = "{wd}/{strain}/.cat_process_{category}_done",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    log:
        "{wd}/{strain}/logs/stage2_process_cat_{category}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/SR_process_categories.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.category} \
            {params.category_path} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_order_categories:
    """
    Order and align category fingerprints across timepoints for each
    pairwise comparison suffix.  Runs after all 13 categories are processed.
    One rule instance per (strain, analysis_suffix).
    """
    input:
        all_cats_done = expand(
            "{{wd}}/{{strain}}/.cat_process_{category}_done",
            category=CATEGORIES,
        ),
    output:
        sentinel = "{wd}/{strain}/.cat_order_{analysis_suffix}_done",
    params:
        scripts_dir   = get_scripts_dir(),
        category_path = get_categories_dir(),
    wildcard_constraints:
        analysis_suffix = "T0vsTSG|T0vsTLG|T0vsTLR",
    log:
        "{wd}/{strain}/logs/stage2_order_cat_{analysis_suffix}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/SR_order_categories.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {params.category_path} \
            {wildcards.analysis_suffix} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_plot_categories:
    """
    Generate per-suffix category comparison SVG plots.
    One rule instance per (strain, analysis_suffix).
    """
    input:
        order_done = "{wd}/{strain}/.cat_order_{analysis_suffix}_done",
    output:
        sentinel = "{wd}/{strain}/.cat_plot_{analysis_suffix}_done",
    params:
        scripts_dir = get_scripts_dir(),
    wildcard_constraints:
        analysis_suffix = "T0vsTSG|T0vsTLG|T0vsTLR",
    log:
        "{wd}/{strain}/logs/stage2_plot_cat_{analysis_suffix}.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/SR_plot_categories.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            {wildcards.analysis_suffix} \
            > {log} 2>&1
        touch {output.sentinel}
        """


rule r_plot_gal_vs_raf:
    """
    Plot galactose vs raffinose coverage comparison.
    Terminal rule for stage 2; produces the .stage2_done sentinel.
    """
    input:
        plots_done = expand(
            "{{wd}}/{{strain}}/.cat_plot_{analysis_suffix}_done",
            analysis_suffix=ANALYSIS_SUFFIXES,
        ),
    output:
        sentinel = "{wd}/{strain}/.stage2_done",
    params:
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/{strain}/logs/stage2_gal_vs_raf.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        Rscript {params.scripts_dir}/R/SR_plot_gal_vs_raf.R \
            {wildcards.wd}/{wildcards.strain}/ \
            {wildcards.strain} \
            > {log} 2>&1

        # Organise outputs into subfolders
        mkdir -p {wildcards.wd}/{wildcards.strain}/Category_Data_TSV_and_Plots
        mv {wildcards.wd}/{wildcards.strain}/*fingerprint* \
           {wildcards.wd}/{wildcards.strain}/Category_Data_TSV_and_Plots/ 2>/dev/null || true
        mv {wildcards.wd}/{wildcards.strain}/*.svg \
           {wildcards.wd}/{wildcards.strain}/Category_Data_TSV_and_Plots/ 2>/dev/null || true

        touch {output.sentinel}
        """
