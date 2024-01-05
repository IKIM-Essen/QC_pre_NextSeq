rule snakemake_report:
    input:
        # 1. Quality control
        "results/{date}/report/qc/multiqc.html",
        "results/{date}/report/plots/human_contamination.html",
        # 2. species diversity
        "results/{date}/report/diversity/diversity_summary/",
        expand(
            "results/{{date}}/report/plots/abundance_{level}.html",
            level=get_bacterial_levels(),
        ),
    output:
        "results/{date}/report/{date}_report.zip",
    log:
        "logs/{date}/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --nolock --report {output} "
        "> {log} 2>&1"
