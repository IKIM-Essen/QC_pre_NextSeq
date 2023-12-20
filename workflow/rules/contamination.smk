rule download_human_ref:
    output:
        get_human_ref(),
    params:
        download=config["human-ref"],
        folder=get_resource_path(),
    log:
        "logs/human_ref_download.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.folder} && "
        "cd {params.folder} && "
        "wget {params.download}) > {log} 2>&1"


rule minimap2_bam_sorted:
    input:
        target=get_human_ref(),  # can be either genome index or genome fasta
        query=get_trimmed_fastqs,
    output:
        "results/{date}/contamination/{sample}.sorted.bam",
    log:
        "logs/{date}/contamination/mapping_{sample}.log",
    params:
        extra="-x map-pb",  # optional
        sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 4
    wrapper:
        "v3.3.1/bio/minimap2/aligner"


rule samtools_index:
    input:
        "results/{date}/contamination/{sample}.sorted.bam",
    output:
        temp("results/{date}/contamination/{sample}.sorted.bam.bai"),
    log:
        "logs/{date}/contamination/indexing_{sample}.log",
    params:
        extra="",
    threads: 4
    wrapper:
        "v3.3.1/bio/samtools/index"


rule host_stats:
    input:
        bam="results/{date}/contamination/{sample}.sorted.bam",
    output:
        "results/{date}/contamination/stats_{sample}.txt",
    params:
        extra="",
    log:
        "logs/{date}/contamination/stats_{sample}.log",
    wrapper:
        "v3.3.1/bio/samtools/stats"


"""

rule map_to_host:
    input:
        fastqs=get_trimmed_fastqs,
    output:
        temp("results/{date}/host_filtering/alignments/{sample}.bam"),
    params:
        ref=config["host_filtering"]["ref_genome"],
    threads: 10
    log:
        "logs/{date}/host_filtering/map_to_host_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "(minimap2 -a -xsr -t {threads} {params.ref} {input.fastqs} | "
        "samtools view -bh | "
        "samtools sort --threads {threads} -o {output}) > {log} 2>&1"

rule index_host_alignment:
    input:
        "results/{date}/host_filtering/{sample}.bam",
    output:
        temp("results/{date}/host_filtering/{sample}.bam.bai"),
    threads: 3
    log:
        "logs/{date}/host_filtering/index_host_alignment_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"

rule filter_host:
    input:
        bam="results/{date}/host_filtering/alignments/{sample}.bam",
        bai="results/{date}/host_filtering/alignments/{sample}.bam.bai",
    output:
        non_host=expand("results/{{date}}/host_filtering/non_host/{{sample}}_{reads}.fastq.gz",reads=["R1","R2"]),
    threads: 3
    log:
        "logs/{date}/host_filtering/filter_host_{sample}.log", 
    conda:
        "../envs/minimap2.yaml"
    shell:
        "(samtools fastq -F 3584 -f 77 {input.bam} | gzip -c > {output.non_host[0]} && "
        "samtools fastq -F 3584 -f 141 {input.bam} | gzip -c > {output.non_host[1]}) > {log} 2>&1"


rule host_filtering_summary:
    input:
        csv="results/{date}/report/kraken2_summary.csv",
        jsons=expand(
            "results/{{date}}/trimmed/fastp/{sample}.fastp.json",
            sample=get_samples(),
        ),
    output:
        csv="results/{date}/report/host_filtering_summary.csv",
    params:
        host_name=config["host_filtering"]["host_name"]
    log:
        "logs/{date}/host_filtering/summary.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/host_filtering_summary.py"

"""
