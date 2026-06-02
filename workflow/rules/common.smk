import os


configfile: "config/config.yaml"


def get_resource_path():
    return config["resources"]


def get_run_date():
    return config["run-date"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards):
    file_r1 = pep.sample_table.loc[wildcards.sample]["fq1"]
    file_r2 = pep.sample_table.loc[wildcards.sample]["fq2"]
    return (
        file_r1,
        file_r2,
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
    if config["human-ref"]["use-local"]:
        local_ref = config["human-ref"]["local-path"]
    else:
        path = config["human-ref"]["download-path"]
        local_ref = "{}{}".format(get_resource_path(), path.split("/")[-1])
    return local_ref


def get_kaiju_fmi_file():
    return config["kaiju-db"]["fmi-file"]


def get_kaiju_nodes_file():
    return config["kaiju-db"]["nodes-dmp"]


def get_kaiju_names_file():
    return config["kaiju-db"]["names-dmp"]


def get_tax_levels():
    return ["genus", "domain"]
