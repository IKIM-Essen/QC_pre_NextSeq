if config["human-ref"]["use-local"]:
    rule copy_local_human_ref:
        output:
            fasta=temp(get_human_ref()),
        params:
            local=config["human-ref"]["local-path"],
            folder=lambda wildcards, output: Path(output.fasta).parent,  #get_resource_path(),
        log:
            "logs/human_ref_local_copy.log",
        group:
            "refGenome_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder} && "
            "cp {params.local} {output.fasta}) > {log} 2>&1"

else:
    rule download_human_ref:
        output:
            fasta=get_human_ref(),
        params:
            download=config["human-ref"]["download-path"],
            folder=lambda wildcards, output: Path(output.fasta).parent,  #get_resource_path(),
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
