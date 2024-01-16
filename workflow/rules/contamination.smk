rule download_human_ref:
    output:
        hfile=get_human_ref(),
    params:
        download=config["human-ref"],
        folder=lambda wildcards, output: Path(output.hfile).parent,  #get_resource_path(),
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
        target=get_human_ref(),
        query=get_trimmed_fastqs,
    output:
        temp("results/{date}/contamination/{sample}.sorted.bam"),
    log:
        "logs/{date}/contamination/mapping/{sample}.log",
    params:
        extra="-x map-pb",
        sorting="coordinate",
        sort_extra="",
    threads: 12
    wrapper:
        "v3.3.3/bio/minimap2/aligner"


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
        "v3.3.3/bio/samtools/stats"
