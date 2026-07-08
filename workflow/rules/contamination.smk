if not config["human-ref"]["use-local"]:

    rule download_human_ref:
        output:
            fasta=get_human_ref(),
        params:
            download=config["human-ref"]["download-path"],
            folder=lambda wildcards, output: Path(output.fasta).parent,
        log:
            "logs/human_ref_download.log",
        group:
            "refGenome_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder} && "
            "cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


rule minimap2_bam_sorted:
    input:
        target=get_human_ref(),
        query=get_trimmed_fastqs,
    output:
        temp("results/{date}/contamination/{sample}.sorted.bam"),
    log:
        "logs/{date}/contamination/mapping/{sample}.log",
    group:
        "refGenome_depended"
    params:
        extra="-x map-sr",
        sorting="coordinate",
        sort_extra="",
    threads: 12
    wrapper:
        "v9.4.2/bio/minimap2/aligner"


rule host_stats:
    input:
        bam=rules.minimap2_bam_sorted.output,
        #"results/{date}/contamination/{sample}.sorted.bam",
    output:
        temp("results/{date}/contamination/{sample}_stats.txt"),
    params:
        extra="",
    log:
        "logs/{date}/contamination/stats/{sample}.log",
    wrapper:
        "v9.4.2/bio/samtools/stats"
