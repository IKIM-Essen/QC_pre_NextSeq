rule qc_diversity_summary:
    input:
        jsons=expand(
            "results/{{date}}/qc/fastp/{sample}.fastp.json",
            sample=get_samples(),
        ),
        stats=expand(
            "results/{{date}}/contamination/{sample}_stats.txt", sample=get_samples()
        ),
        bracken="results/{date}/report/bracken/merged.bracken_domain.txt",
    output:
        summary_csv=ensure(
            "results/{date}/report/filtering_summary.csv", non_empty=True
        ),
        human_cont_html=report(
            "results/{date}/report/plots/human_contamination.html",
            category="5. Human contamination plot",
        ),
        read_summary_html=report(
            "results/{date}/report/plots/read_summary.html",
            category="2. Number of reads plot",
        ),
        domain_abd_html=report(
            "results/{date}/report/plots/domain_abundance.html",
            category="4. Domain level abundance plot",
        ),
    log:
        "logs/{date}/summary_and_plots.log",
    threads: 4
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_plots.py"


rule summary2report:
    input:
        "results/{date}/report/filtering_summary.csv",
    output:
        report(
            directory("results/{date}/report/filtering_summary/"),
            htmlindex="index.html",
            category="1. Summary table",
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="summary",
        header=" ",
        pattern=config["tablular-config"],
    log:
        "logs/{date}/summary_to_html.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --pin-until {params.pin_until} {output} && "
        "(sed -i '{params.pattern} {params.header}</a>' "
        "{output}/indexes/index1.html && "
        "sed -i 's/report.xlsx/{params.name}_report.xlsx/g' {output}/indexes/index1.html) && "
        "mv {output}/report.xlsx {output}/{params.name}_report.xlsx && "
        "cp {params.styles}* {output}/css/ > {log} 2>&1"


if not config["testing"]:

    rule snakemake_report:
        input:
            "results/{date}/report/filtering_summary/",
            "results/{date}/report/qc/multiqc.html",
            rules.qc_diversity_summary.output.read_summary_html,
            rules.qc_diversity_summary.output.human_cont_html,
            rules.qc_diversity_summary.output.domain_abd_html,
        output:
            "results/{date}/report/{date}_report.zip",
        log:
            "logs/{date}/snakemake-report.log",
        conda:
            "../envs/snakemake.yaml"
        shell:
            "snakemake --nolock --report {output} --profile '' > {log} 2>&1"
