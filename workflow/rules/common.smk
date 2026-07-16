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


def _cfg_bool(key, default):
    """A config flag as a bool, tolerating a real YAML bool or a string
    (snakemake --config passes strings)."""
    val = config.get(key, default)
    if isinstance(val, str):
        return val.strip().lower() in ("true", "1", "yes", "on")
    return bool(val)


def emit_preprocessed():
    """Whether to PRODUCE + deposit the preprocessed reads at all. run_qc.sh sets
    it False for 16S (which keeps its own QIIME2 QC — primer trimming + DADA2 — so
    fastp-trimmed reads aren't consumable there; only QC metrics are wanted) and
    True for isolate/metagenome. Default True."""
    return _cfg_bool("emit-preprocessed", True)


def remove_human():
    """Whether to REMOVE human reads from the deposited preprocessed reads (not
    just MEASURE %human — that always happens). Type-dependent in practice:
    run_qc.sh sets it True for isolate/metagenome and False for 16S (amplicon,
    where host filtering isn't meaningful). Accepts a real YAML bool or a string
    ('--config remove-human=False' passes a string), defaulting to True."""
    return _cfg_bool("remove-human", True)


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
    # In QC mode kaiju is skipped entirely, so no taxonomic levels are
    # requested -> every expand(..., level=get_tax_levels()) becomes empty and
    # the kaiju targets drop out of `rule all`. In diversity mode kaiju runs.
    if config.get("mode", "qc") == "qc":
        return []
    return ["genus", "domain"]
