rule snakemake_report:
    input:
        # 1. Quality control
        "results/{date}/output/multiqc.html",
        "results/{date}/output/plot_contamination.html",
    output:
        "results/{date}/output/report.zip",  #html",
    log:
        "logs/{date}/snakemake-report.log",
    conda:
        "../envs/snakemake.yaml"
    shell:
        "snakemake --nolock --report {output} "
        "> {log} 2>&1"
