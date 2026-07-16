rule qc_diversity_summary:
    input:
        jsons=expand(
            "results/{{date}}/qc/fastp/{sample}.fastp.json",
            sample=get_samples(),
        ),
        stats=expand(
            "results/{{date}}/contamination/{sample}_stats.txt", sample=get_samples()
        ),
        nonhuman=expand(
            "results/{{date}}/contamination/{sample}_nonhuman_pairs.txt",
            sample=get_samples(),
        ),
        kaiju=expand(
            "results/{{date}}/report/kaiju/merged.kaiju_{level}.tsv",
            level=get_tax_levels(),
        ),
    output:
        summary_csv="results/{date}/report/filtering_summary.csv",
        # machine-readable per-sample metrics for the sample registry (raw
        # numbers); ingested by dispatcher.py --reconcile into qc_metrics.
        qc_metrics_tsv="results/{date}/report/qc_metrics.tsv",
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
        genus_abd_html=report(
            "results/{date}/report/plots/genus_abundance.html",
            category="5. Genus level abundance plot",
        ),
        genus_top10_csv="results/{date}/report/kaiju/genus_top10.csv",
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


rule genus_top10_report:
    input:
        "results/{date}/report/kaiju/genus_top10.csv",
    output:
        report(
            directory("results/{date}/report/genus_top10/"),
            htmlindex="index.html",
            category="6. Genus top 10 table",
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="genus_top10",
        header=" ",
        pattern=config["tablular-config"],
    log:
        "logs/{date}/genus_top10_to_html.log",
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
            # kaiju-based sections only exist meaningfully in diversity mode;
            # in QC mode create_plots writes placeholders to satisfy the DAG,
            # but we must NOT bundle them into the report (they were the "fake"
            # empty domain/genus/top10 chapters). Include them only for diversity.
            *(
                [
                    rules.qc_diversity_summary.output.domain_abd_html,
                    rules.qc_diversity_summary.output.genus_abd_html,
                    rules.genus_top10_report.output,
                ]
                if config.get("mode", "qc") == "diversity"
                else []
            ),
        output:
            "results/{date}/report/{date}_report.zip",
        log:
            "logs/{date}/snakemake-report.log",
        conda:
            "../envs/snakemake.yaml"
        shell:
            "snakemake --nolock --report {output} > {log} 2>&1"
