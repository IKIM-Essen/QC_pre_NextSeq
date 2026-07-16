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
        extra="-x sr",
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


rule nonhuman_pairs:
    # Count the read PAIRS that survive human removal (both mates unmapped). This
    # is exactly the R1 selection of the host-removal step (samtools fastq -f 77),
    # so -c on the same filter == number of retained pairs. Reported as
    # host_removed_read_pairs in qc_metrics. Computed for ALL types (the bam
    # always exists for the %human metric); for 16S no removal actually happens,
    # so it is the informational count of non-human pairs.
    input:
        bam=rules.minimap2_bam_sorted.output,
    output:
        temp("results/{date}/contamination/{sample}_nonhuman_pairs.txt"),
    log:
        "logs/{date}/contamination/nonhuman_pairs/{sample}.log",
    threads: 2
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -c -F 3584 -f 77 {input.bam} > {output} 2> {log}"
