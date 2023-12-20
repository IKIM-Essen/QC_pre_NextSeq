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


def get_trimmed_fastq(wildcards):
    path = get_data_path()
    return (
        "{data}{{date}}/{{sample}}_R1.fastq.gz".format(data=path),
    )


def get_adapters(wildcards):
    return config["adapter-seqs"]


def get_trimmed_fastqs(wildcards):
    return [
        "results/{date}/qc/fastp/{sample}.1.fastq.gz",
        "results/{date}/qc/fastp/{sample}.2.fastq.gz",
    ]


def get_human_ref():
    link = config["human-ref"]
    path = get_resource_path()
    file = link.split("/")[-1]
    return f"{path}{file}"
