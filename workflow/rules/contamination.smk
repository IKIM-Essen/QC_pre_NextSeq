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
        target=get_human_ref(),
        query=get_trimmed_fastqs,
    output:
        temp("results/{date}/contamination/{sample}.sorted.bam"),
    log:
        "logs/{date}/contamination/mapping_{sample}.log",
    params:
        extra="-x map-pb",  
        sorting="coordinate",
        sort_extra="",
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


rule plot_contamination:
    input:
        stats=expand(
            "results/{{date}}/contamination/stats_{sample}.txt", sample=get_samples()
        ),
    output:
        report(
            "results/{date}/output/plot_contamination.html",
            category="Quality control",
        ),
    log:
        "logs/{date}/contamination/plot.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/cont_calc.py"
