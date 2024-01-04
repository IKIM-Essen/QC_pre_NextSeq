rule snakemake_report:
    input:
        # 1. Quality control
        "results/{date}/qc/multiqc.html",
        "results/{date}/report/contamination.html",
        # 2. species diversity
        "results/{date}/report/diversity/",
        expand("results/{{date}}/report/{level}_abundance.html", level=get_bacterial_levels())
    output:
        "results/{date}/report/report.zip",  #html",
    log:
        "logs/{date}/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --nolock --report {output} "
        "> {log} 2>&1"
