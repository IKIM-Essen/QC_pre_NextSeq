from pathlib import Path


def get_kaiju_fmi_file():
    if config["kaiju-db"]["use-local"]:
        return config["kaiju-db"]["local-fmi"]
    return config["kaiju-db"]["fmi"]


def get_kaiju_nodes_file():
    if config["kaiju-db"]["use-local"]:
        return config["kaiju-db"]["local-nodes-dmp"]
    return config["kaiju-db"]["nodes-dmp"]


def get_kaiju_names_file():
    if config["kaiju-db"]["use-local"]:
        return config["kaiju-db"]["local-names-dmp"]
    return config["kaiju-db"]["names-dmp"]


if not config["kaiju-db"]["use-local"]:

    rule download_kaiju_db:
        output:
            fmi=config["kaiju-db"]["fmi"],
            nodes=config["kaiju-db"]["nodes-dmp"],
            names=config["kaiju-db"]["names-dmp"],
        params:
            download=config["kaiju-db"]["download-path"],
            db_folder=lambda wc, output: Path(output.fmi).parent,
        log:
            "logs/kaiju_DB_download.log",
        group:
            "kaijuDB_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p '{params.db_folder}' && "
            "wget -c '{params.download}' -O - | "
            "tar -xzv -C '{params.db_folder}') > {log} 2>&1"


rule kaiju:
    input:
        fmi=get_kaiju_fmi_file(),
        nodes=get_kaiju_nodes_file(),
        names=get_kaiju_names_file(),
        fastqs=get_trimmed_fastqs,
    output:
        kout=temp("results/{date}/diversity/kaiju_outfiles/{sample}.out"),
    threads: 12
    log:
        "logs/{date}/kaiju/run/{sample}.log",
    group:
        "kaijuDB_depended"
    conda:
        "../envs/kaiju_based.yaml"
    shell:
        "kaiju -t {input.nodes} -f {input.fmi} "
        "-i {input.fastqs[0]} -j {input.fastqs[1]} "
        "-z {threads} -o {output.kout} > {log} 2>&1"


rule kaiju_genus:
    input:
        nodes=get_kaiju_nodes_file(),
        names=get_kaiju_names_file(),
        kout=rules.kaiju.output.kout,
    output:
        report=temp("results/{date}/report/kaiju/reports_genus/{sample}.tsv"),
    params:
        rank="genus",
    log:
        "logs/{date}/kaiju/genus/{sample}.log",
    threads: 2
    group:
        "kaijuDB_depended"
    conda:
        "../envs/kaiju_based.yaml"
    shell:
        "kaiju2table -t {input.nodes} -n {input.names} "
        "-r {params.rank} -o {output.report} {input.kout} > {log} 2>&1"


use rule kaiju_genus as kaiju_domain with:
    output:
        report=temp("results/{date}/report/kaiju/reports_domain/{sample}.tsv"),
    params:
        rank="superkingdom",
    log:
        "logs/{date}/kaiju/domain/{sample}.log",


rule merge_kaiju:
    input:
        expand(
            "results/{{date}}/report/kaiju/reports_{{level}}/{sample}.tsv",
            sample=get_samples(),
        ),
    output:
        "results/{date}/report/kaiju/merged.kaiju_{level}.tsv",
    log:
        "logs/{date}/kaiju/merge_{level}.log",
    threads: 1
    conda:
        "../envs/kaiju_based.yaml"
    shell:
        "cat {input} > {output} 2> {log}"
