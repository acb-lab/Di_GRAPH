"""
Stage 5 — HTML Report Generation.

Renders the Di-GRAPH flexdashboard report from ``Di-GRAPH_report.Rmd``
using ``generate_report.R``.  This is the terminal stage of the pipeline;
it waits for every per-strain stage 1–4 sentinel before running.

R script interface:
  generate_report.R   (no arguments — renders Di-GRAPH_report.Rmd from
                       the current working directory and writes
                       Di-GRAPH_report.html in that same directory)

The rule runs in the report files directory so that ``render()`` can locate
the .Rmd file by its basename, then moves the generated HTML to the
experiment working directory.
"""


rule generate_report:
    """
    Render the Di-GRAPH flexdashboard HTML report.
    Terminal rule for the full pipeline.
    """
    input:
        stage1_done = expand(
            "{{wd}}/{strain}/.stage1_done",
            strain=get_sample_names(),
        ),
        stage2_done = expand(
            "{{wd}}/{strain}/.stage2_done",
            strain=get_sample_names(),
        ),
        stage3_done = expand(
            "{{wd}}/{strain}/.stage3_done",
            strain=get_sample_names(),
        ),
        stage4_done = expand(
            "{{wd}}/{strain}/.stage4_done",
            strain=get_sample_names(),
        ),
    output:
        report = "{wd}/Di-GRAPH_report.html",
    params:
        report_dir  = get_report_dir(),
        scripts_dir = get_scripts_dir(),
    log:
        "{wd}/logs/stage5_report.log"
    conda:
        "../envs/digraph.yaml"
    shell:
        """
        cd {params.report_dir}
        Rscript generate_report.R > {log} 2>&1
        mv Di-GRAPH_report.html {wildcards.wd}/Di-GRAPH_report.html
        """
