rule download_kraken_db:
    output:
        hfile=get_kraken_db_file(),
    params:
        download=config["kraken"]["download-path"],
        db_folder=lambda wildcards, output: Path(output.hfile).parent,
    log:
        "logs/kraken2_DB_download.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.db_folder} && "
        "wget -c {params.download} -O - | "
        "tar -zxv -C {params.db_folder}) > {log} 2>&1"


rule kraken2:
    input:
        hfile=get_kraken_db_file(),
        fastqs=get_trimmed_fastqs,
    output:
        report="results/{date}/diversity/kraken_reports/{sample}_report.tsv",
        outfile=temp("results/{date}/diversity/kraken_outfiles/{sample}_outfile.tsv"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
    threads: 32
    log:
        "logs/{date}/diversity/kraken2_run/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --quick --paired "
        "--output {output.outfile} --report {output.report} "
        "--gzip-compressed {input.fastqs} > {log} 2>&1"


rule kraken_summary:
    input:
        reports=expand(
            "results/{{date}}/diversity/kraken_reports/{sample}_report.tsv",
            sample=get_samples(),
        ),
        jsons=expand(
            "results/{{date}}/qc/fastp/{sample}.fastp.json",
            sample=get_samples(),
        ),
    output:
        csv="results/{date}/diversity/diversity_summary.csv",
    log:
        "logs/{date}/diversity/summary.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/kraken_summary.py"


rule kraken2_report:
    input:
        "results/{date}/diversity/diversity_summary.csv",
    output:
        report(
            directory("results/{date}/report/diversity/"),
            htmlindex="index.html",
            category="2. Species diversity",
            labels={"level":"all"},
        ),
    params:
        pin_until="sample",
        styles="resources/report/tables/",
        name="diversity_summary",
        header="Kraken2 diversity summary",
    log:
        "logs/{date}/diversity/summary_to_html.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --pin-until {params.pin_until} {output} && "
        "(sed -i '/>github<\/a>/a \\\\t\\t\\t</li>\\n\\t\\t\\t<li class=\"nav-item\">"
        "\\n\\t\\t\\t\\t<a class=\"nav-link\" href=\"#\">{params.header}</a>' "
        "{output}/indexes/index1.html && "
        "sed -i 's/report.xlsx/{params.name}_report.xlsx/g' {output}/indexes/index1.html) && "
        "mv {output}/report.xlsx {output}/{params.name}_report.xlsx && "
        "cp {params.styles}* {output}/css/ > {log} 2>&1"


rule bracken_genus:
    input:
        hfile=get_kraken_db_file(),
        kreport=get_kraken_report,
    output:
        breport=temp("results/{date}/report/bracken/reports_genus/{sample}.breport"),
        bfile=temp("results/{date}/report/bracken/files_genus/{sample}.bracken"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "G",
    log:
        "logs/{date}/diversity/bracken/genus/{sample}.log",
    resources:
        mem_mb=100,
    threads: 1
    conda:
        "../envs/bracken.yaml"
    shell:
        "bracken -d {params.db} -i {input.kreport} -l {params.level} -o {output.bfile} -w {output.breport} > {log} 2>&1"


use rule bracken_genus as bracken_family with:
    input:
        hfile=get_kraken_db_file(),
        kreport=get_kraken_report,
    output:
        breport=temp("results/{date}/report/bracken/reports_family/{sample}.breport"),
        bfile=temp("results/{date}/report/bracken/files_family/{sample}.bracken"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "F",
    log:
        "logs/{date}/diversity/bracken/family/{sample}.log",


use rule bracken_genus as bracken_phylum with:
    input:
        hfile=get_kraken_db_file(),
        kreport=get_kraken_report,
    output:
        breport=temp("results/{date}/report/bracken/reports_phylum/{sample}.breport"),
        bfile=temp("results/{date}/report/bracken/files_phylum/{sample}.bracken"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "P",
    log:
        "logs/{date}/diversity/bracken/phylum/{sample}.log",


use rule bracken_genus as bracken_class with:
    input:
        hfile=get_kraken_db_file(),
        kreport=get_kraken_report,
    output:
        breport=temp("results/{date}/report/bracken/reports_class/{sample}.breport"),
        bfile=temp("results/{date}/report/bracken/files_class/{sample}.bracken"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level= "C",
    log:
        "logs/{date}/diversity/bracken/class/{sample}.log",


rule merge_bracken:
    input:
        expand(
            "results/{{date}}/report/bracken/files_{{level}}/{sample}.bracken",
            sample=get_samples(),
        ),
    output:
        "results/{date}/report/bracken/merged.bracken_{level}.txt",
    log:
        "logs/{date}/diversity/bracken/merge_{level}.log",
    resources:
        mem_mb=100,
    threads: 1
    params:
        threads=1,
    conda:
        "../envs/bracken.yaml"
    shell:
        "(python $CONDA_PREFIX/bin/combine_bracken_outputs.py --files {input} --output {output}) > {log} 2>&1"


'''
rule create_bracken_plot:
    input:
        "results/{date}/report/bracken/merged.bracken_{level}.txt",
    output:
        report(
            "results/{date}/report/bracken_{level}_plot.png",
            category="2. Species diversity",
            labels={"level":"{level}"},
        ),
    params:
        threshold=0.001,
    log:
        "logs/{date}/diversity/bracken/plot_{level}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/brackenplot.py"
'''