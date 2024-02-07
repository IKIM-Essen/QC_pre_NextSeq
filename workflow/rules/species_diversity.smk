if config["kraken-db"]["use-local"]:

    rule copy_local_kraken_db:
        output:
            hfile=get_kraken_db_file(),
        params:
            local=config["kraken-db"]["local-path"],
            db_folder=lambda wildcards, output: Path(output.hfile).parent,
            resource_folder=lambda wildcards, output: Path(output.hfile).parent.parent,
            filename=get_kraken_db_tar(),
        log:
            "logs/kraken2_DB_local_copy.log",
        group:
            "krakenDB_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.db_folder}/ && "
            "cp {params.local} {params.resource_folder}/ && "
            "tar fzxv {params.resource_folder}/{params.filename} -C {params.db_folder}/ && "
            "rm {params.resource_folder}/{params.filename}) > {log} 2>&1"

else:

    rule download_kraken_db:
        output:
            hfile=get_kraken_db_file(),
        params:
            download=get_kraken_db_url(),
            db_folder=lambda wildcards, output: Path(output.hfile).parent,
        log:
            "logs/kraken2_DB_download.log",
        group:
            "krakenDB_depended"
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
        report=temp("results/{date}/diversity/kraken_reports/{sample}_report.tsv"),
        outfile=temp("results/{date}/diversity/kraken_outfiles/{sample}_outfile.tsv"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
    threads: 32
    log:
        "logs/{date}/kraken2_run/{sample}.log",
    group:
        "krakenDB_depended"
    conda:
        "../envs/kraken_based.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --quick --paired "
        "--output {output.outfile} --report {output.report} "
        "--gzip-compressed {input.fastqs} > {log} 2>&1"


rule bracken_genus:
    input:
        hfile=get_kraken_db_file(),
        kreport=get_kraken_report,
    output:
        breport=temp("results/{date}/report/bracken/reports_genus/{sample}.breport"),
        bfile=temp("results/{date}/report/bracken/files_genus/{sample}.bracken"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level="G",
    log:
        "logs/{date}/bracken/genus/{sample}.log",
    resources:
        mem_mb=100,
    threads: 2
    group:
        "krakenDB_depended"
    conda:
        "../envs/kraken_based.yaml"
    shell:
        "bracken -d {params.db} -i {input.kreport} -l {params.level} -o {output.bfile} -w {output.breport} > {log} 2>&1"


use rule bracken_genus as bracken_domain with:
    input:
        hfile=get_kraken_db_file(),
        kreport=get_kraken_report,
    output:
        breport=temp("results/{date}/report/bracken/reports_domain/{sample}.breport"),
        bfile=temp("results/{date}/report/bracken/files_domain/{sample}.bracken"),
    params:
        db=lambda wildcards, input: Path(input.hfile).parent,
        level="D",
    log:
        "logs/{date}/bracken/domain/{sample}.log",
    group:
        "krakenDB_depended"


rule merge_bracken:
    input:
        expand(
            "results/{{date}}/report/bracken/files_{{level}}/{sample}.bracken",
            sample=get_samples(),
        ),
    output:
        "results/{date}/report/bracken/merged.bracken_{level}.txt",
    log:
        "logs/{date}/bracken/merge_{level}.log",
    resources:
        mem_mb=100,
    threads: 2
    params:
        threads=1,
    conda:
        "../envs/kraken_based.yaml"
    shell:
        "(python $CONDA_PREFIX/bin/combine_bracken_outputs.py --files {input} --output {output}) > {log} 2>&1"
