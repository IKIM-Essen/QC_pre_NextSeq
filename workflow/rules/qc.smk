RAW_DATA_PATH = get_data_path()


rule local_fastqs:
    input:
        get_fastqs,
    output:
        raw1=temp(f"{RAW_DATA_PATH}{date}/{sample}_R1.fastq.gz"),
        raw2=temp(f"{RAW_DATA_PATH}{date}/{sample}_R2.fastq.gz"),
    params:
        outdir=lambda wildcards, output: Path(output.raw1).parent,
    log:
        "logs/{date}/copy_data/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.outdir} && "
        "cp -v -t {params.outdir} {input}) > {log} 2>&1"


rule fastp:
    input:
        sample=get_local_fastqs,
    output:
        trimmed=temp(
            [
                "results/{date}/qc/fastp/{sample}.1.fastq.gz",
                "results/{date}/qc/fastp/{sample}.2.fastq.gz",
            ]
        ),
        html=temp("results/{date}/qc/fastp/{sample}.html"),
        json=temp("results/{date}/qc/fastp/{sample}.fastp.json"),
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {phred} --length_required {minlen}".format(
            phred=(config["quality-criteria"]["min-PHRED"]),
            minlen=(config["quality-criteria"]["min-length-reads"]),
        ),
    log:
        "logs/{date}/qc/fastp/{sample}.log",
    threads: 2
    wrapper:
        "v3.3.1/bio/fastp"


rule fastqc:
    input:
        get_trimmed_fastq,
    output:
        html=temp("results/{date}/qc/fastqc/{sample}.html"),
        zip=temp("results/{date}/qc/fastqc/{sample}_fastqc.zip"),
    log:
        "logs/{date}/qc/fastqc/{sample}.log",
    threads: 4
    resources:
        mem_mb=1024,
    wrapper:
        "v3.3.1/bio/fastqc"


rule multiqc:
    input:
        expand(
            [
                "results/{{date}}/qc/fastqc/{sample}_fastqc.zip",
                "results/{{date}}/qc/fastp/{sample}.fastp.json",
            ],
            sample=get_samples(),
        ),
    output:
        report(
            "results/{date}/report/qc/multiqc.html",
            category="1. Quality control",
        ),
        "results/{date}/report/qc/multiqc_data.zip",
    params:
        extra=(
            "--zip-data-dir "
            "--config config/multiqc_config.yaml "
            "--title 'Quality control for MiSeq run from {date}'"
        ),
        use_input_files_only=True,
    log:
        "logs/{date}/qc/multiqc.log",
    wrapper:
        "v3.3.1/bio/multiqc"
