rule fastp:
    input:
        sample=get_fastqs,
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
        "v9.8.0/bio/fastp"


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
        "v7.6.0/bio/fastqc"


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
            category="3. MultiQC report",
        ),
        #labels={"": "3. MultiQC report"}
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
        "v9.8.0/bio/multiqc"
