if config["kraken-db"]["use-local"]:

    rule copy_local_kraken_db:
        output:
            hfile=get_kraken_db_file(),
        params:
            local=config["kraken-db"]["local-path"],
            db_folder=lambda wildcards, output: Path(output.hfile).parent,
            resource_folder=lambda wildcards, output: Path(output.hfile).parent.parent,
            filename=get_kraken_db_tar(),
        log:
            "logs/kraken2_DB_local_copy.log",
        group:
            "krakenDB_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.db_folder}/ && "
            "cp {params.local} {params.resource_folder}/ && "
            "tar fzxv {params.resource_folder}/{params.filename} -C {params.db_folder}/ && "
            "rm {params.resource_folder}/{params.filename}) > {log} 2>&1"

else:

    rule download_kraken_db:
        output:
            hfile=get_kraken_db_file(),
        params:
            download=get_kraken_db_url(),
            db_folder=lambda wildcards, output: Path(output.hfile).parent,
        log:
            "logs/kraken2_DB_download.log",
        group:
            "krakenDB_depended"
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.db_folder} && "
            "wget -c {params.download} -O - | "
            "tar -zxv -C {params.db_folder}) > {log} 2>&1"


#rule kraken2:
#    input:
 #       hfile=get_kraken_db_file(),
  #      fastqs=get_trimmed_fastqs,
   # output:
    #    report=temp("results/{date}/diversity/kraken_reports/{sample}_report.tsv"),
     #   outfile=temp("results/{date}/diversity/kraken_outfiles/{sample}_outfile.tsv"),
    #params:
     #   db=lambda wildcards, input: Path(input.hfile).parent,
    #threads: 32
#   #log:
     #   "logs/{date}/kraken2_run/{sample}.log",
    #group:
     #   "krakenDB_depended"
    #conda:
     #   "../envs/kraken_based.yaml"
    #shell:
    #    "kraken2 --db {params.db} --threads {threads} --quick --paired "
     #   "--output {output.outfile} --report {output.report} "
      #  "--gzip-compressed {input.fastqs} > {log} 2>&1"


# ======================================
# 1. Sketch pro Sample
# ======================================
rule sourmash_sketch:
    input:
        r1="results/{date}/qc/fastp/{sample}.1.fastq.gz",
        r2="results/{date}/qc/fastp/{sample}.2.fastq.gz",
    output:
        sig="results/{date}/report/sourmash/{sample}.sig"
    conda:
        "../envs/sourmash.yaml"
    shell:
        """
        # Erst beide sketchen
        sourmash sketch dna \
          -p k=31,scaled=1000 \
          {input.r1} {input.r2} \
          -o {output.sig}.temp
        
        # Dann mergen zu einer Signatur
        sourmash sig merge \
          {output.sig}.temp \
          -o {output.sig} \
          --name {wildcards.sample}
        
        rm {output.sig}.temp
        """



rule sourmash_sketch_human_db:
    input:
        fasta="resources/GCA_000001405.29_GRCh38.p14_genomic.fna.gz",
    output:
        sig="resources/sourmash_db/human_genome.sig",
    conda:
        "../envs/sourmash.yaml"
    shell:
        "sourmash sketch dna -p k=31,scaled=1000 {input.fasta} -o {output.sig}"

# ======================================
# 3. Gather pro Sample
# ======================================
rule sourmash_gather_sample:
    input:
        sig="results/{date}/report/sourmash/{sample}.sig",
        dbs=[
            "resources/sourmash_db/gtdb-reps-rs226-k31.dna.zip",
            "resources/sourmash_db/ncbi-viruses-2025.01.dna.k=31.sig.zip",
            "resources/sourmash_db/human_genome.sig"
        ]
    output:
        tsv="results/{date}/report/sourmash/{sample}_gather.tsv",
        csv="results/{date}/report/sourmash/{sample}_gather.csv"
    log:
        "logs/{date}/sourmash/gather/{sample}.log"
    conda:
        "../envs/sourmash.yaml"
    shell:
        """
        # WICHTIG: KEIN --rna flag mehr!
        sourmash gather \
            --ksize 31 \
            --threshold-bp 0 \
            {input.sig} \
            {input.dbs} \
            -o {output.csv} \
            > {log} 2>&1 || true
        
        # Erstelle TSV auch wenn keine Matches
        if [ -f {output.csv} ] && [ -s {output.csv} ]; then
            sed 's/,/\t/g' {output.csv} > {output.tsv}
        else
            # Leere Dateien mit Header
            echo -e "intersect_bp\tf_orig_query\tf_match\tf_unique_to_query\tf_unique_weighted\taverage_abund\tmedian_abund\tstd_abund\tname\tfilename\tmd5\tf_match_orig\tunique_intersect_bp\tgather_result_rank\tremaining_bp\tquery_filename\tquery_name\tquery_md5\tquery_bp\tksize\tmoltype\tscaled\tquery_n_hashes\tsum_weighted_found\ttotal_weighted_hashes" > {output.tsv}
            echo "intersect_bp,f_orig_query,f_match,f_unique_to_query,f_unique_weighted,average_abund,median_abund,std_abund,name,filename,md5,f_match_orig,unique_intersect_bp,gather_result_rank,remaining_bp,query_filename,query_name,query_md5,query_bp,ksize,moltype,scaled,query_n_hashes,sum_weighted_found,total_weighted_hashes" > {output.csv}
        fi
        """


rule sourmash_domain_qc:
    input:
        gather=expand(
            "results/{{date}}/report/sourmash/{sample}_gather.tsv",
            sample=get_samples()
        )
    output:
        "results/{date}/report/sourmash/domain_qc.tsv"
    log:
        "logs/{date}/sourmash/domain_qc.log"
    conda:
        "../envs/sourmash.yaml"
    script:
        "../scripts/sourmash_domain_qc.py"

