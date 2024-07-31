import os

ILLUMINA = "illumina"
ONT = "ont"

configfile: "config/config.yaml"


def get_data_path():
    return config["data-handling"]["data"]


def get_resource_path():
    return config["data-handling"]["resources"]


def get_run_date():
    return config["run-date"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_technology(wildcards, sample=None):
    if sample is None:
        sample = wildcards.sample

    return pep.sample_table.loc[sample]["technology"]


def is_ont(wildcards, sample=None):
    if sample is None:
        return get_technology(wildcards) == ONT
    return get_technology(None, sample) == ONT


def is_illumina(wildcards, sample=None):
    if sample is None:
        return get_technology(wildcards) == ILLUMINA
    return get_technology(None, sample) == ILLUMINA


def get_fastqs(wildcards):
    if is_illumina(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]
    elif is_ont(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1"]]


def get_adapters(wildcards):
    return config["adapter-seqs"]


def get_trimmed_fastqs(wildcards):
    if is_illumina(wildcards):
        return [
            "results/{date}/qc/fastp-pe/{sample}.1.fastq.gz",
            "results/{date}/qc/fastp-pe/{sample}.2.fastq.gz",
        ]
    elif is_ont(wildcards):
        return ["results/{date}/qc/fastp-se/{sample}.fastq.gz"]


def get_fastp_results(wildcards):
    """Returns paths of files to aggregate the fastp results for the multiqc rule."""
    # fastp is only used on Illumina and Ion Torrent data
    files = []
    samples = get_samples()
    for sample in samples:
        if is_illumina(None, sample):
            files.append("results/{{date}}/qc/fastp-pe/{sample}.fastp.json".format(sample=sample))
        elif is_ont(None, sample):
            files.append("results/{{date}}/qc/fastp-se/{sample}.fastp.json".format(sample=sample))
    return files


def get_trimmed_fastq(wildcards):
    fastqs = get_trimmed_fastqs(wildcards)
    return fastqs[0]


def get_human_ref():
    if config["human-ref"]["use-local"]:
        path = config["human-ref"]["local-path"]
    else:
        path = config["human-ref"]["download-path"]
    local_ref = "{}{}".format(get_resource_path(), path.split("/")[-1])
    return local_ref


def get_kraken_db_url():
    return config["kraken-db"]["download-path"]


def get_kraken_db_file():
    if config["kraken-db"]["use-local"]:
        path = config["kraken-db"]["local-path"]
    else:
        path = get_kraken_db_url()
    db_name = (Path(path).name).rsplit("_", 1)[0]
    file = "{}{}/hash.k2d".format(get_resource_path(), db_name)
    return file


def get_kraken_db_tar():
    return Path(config["kraken-db"]["local-path"]).name


def get_kraken_report(wildcards):
    return "results/{date}/diversity/kraken_reports/{sample}_report.tsv"


def get_tax_levels():
    return ["genus", "domain"]
