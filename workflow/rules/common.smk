import os


configfile: "config/config.yaml"


def get_data_path():
    return config["data-handling"]["data"]


def get_resource_path():
    return config["data-handling"]["resources"]


def get_run_date():
    return config["run-date"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards):
    return (
        pep.sample_table.loc[wildcards.sample]["fq1"],
        pep.sample_table.loc[wildcards.sample]["fq2"],
    )


def get_local_fastqs(wildcards):
    path = get_data_path()
    return (
        "{data}{{date}}/{{sample}}_R1.fastq.gz".format(data=path),
        "{data}{{date}}/{{sample}}_R2.fastq.gz".format(data=path),
    )


def get_adapters(wildcards):
    return config["adapter-seqs"]


def get_trimmed_fastqs(wildcards):
    return [
        "results/{date}/qc/fastp/{sample}.1.fastq.gz",
        "results/{date}/qc/fastp/{sample}.2.fastq.gz",
    ]


def get_trimmed_fastq(wildcards):
    fastqs = get_trimmed_fastqs(wildcards)
    return fastqs[0]


def get_human_ref():
    file = config["human-ref"].split("/")[-1]
    local_ref = "{}{}".format(get_resource_path(), file)
    return local_ref


def get_kraken_db_file():
    file = "{}{}/hash.k2d".format(get_resource_path(), config["kraken"]["db-name"])
    return file


def get_kraken_report(wildcards):
    return "results/{date}/diversity/kraken_reports/{sample}_report.tsv"


def get_bacterial_levels():
    return ["genus", "family", "class", "phylum"]
