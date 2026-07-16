# =============================================================================
# preprocessed reads — the reads we hand back to the type pipelines so they can
# skip their own fastp/host-filtering and consume already-preprocessed reads.
#
# Two variants, selected by the `remove-human` config flag:
#   remove-human: True  -> fastp-trimmed AND human-removed. The human alignment
#                          (minimap2_bam_sorted) already exists for the %human
#                          contamination metric; we pull the both-mates-unmapped
#                          (= non-human) pairs back out of it with samtools fastq
#                          — the exact host-removal ResMAG does. This is the
#                          default.
#   remove-human: False -> fastp-trimmed only (human is still MEASURED for the
#                          qc_metrics %human, just not removed from the reads).
#
# Either way the output lands at a stable per-sample path. run_qc.sh copies it
# into <seqdata_base>/<type>/fastq/preprocessed/ (the canonical location the type
# pipelines read from) and points the .qc.complete marker at it via a
# `preprocessed=` field; dispatcher.py --reconcile records it as
# result_artifact(kind='preprocessed'). The type pipelines actually consuming
# these reads (skipping their own fastp) is a separate, later step.
# =============================================================================


if remove_human():

    rule preprocessed_reads:
        input:
            # coordinate-sorted alignment of the fastp-trimmed reads against the
            # human reference, already produced for the contamination metric.
            bam=rules.minimap2_bam_sorted.output,
        output:
            r1="results/{date}/preprocessed/{sample}.1.fastq.gz",
            r2="results/{date}/preprocessed/{sample}.2.fastq.gz",
        log:
            "logs/{date}/preprocessed/{sample}.log",
        threads: 8
        conda:
            "../envs/samtools.yaml"
        shell:
            # -F 3584 drops secondary/supplementary/dup/qcfail; -f 77 / -f 141
            # keep the first / second mate of pairs where BOTH mates are unmapped
            # (i.e. non-human). Piping through gzip is samtools-version-agnostic.
            "(samtools fastq --threads {threads} -F 3584 -f 77 {input.bam} | gzip > {output.r1} && "
            "samtools fastq --threads {threads} -F 3584 -f 141 {input.bam} | gzip > {output.r2}) > {log} 2>&1"

else:

    rule preprocessed_reads:
        input:
            r1="results/{date}/qc/fastp/{sample}.1.fastq.gz",
            r2="results/{date}/qc/fastp/{sample}.2.fastq.gz",
        output:
            r1="results/{date}/preprocessed/{sample}.1.fastq.gz",
            r2="results/{date}/preprocessed/{sample}.2.fastq.gz",
        log:
            "logs/{date}/preprocessed/{sample}.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cp {input.r1} {output.r1} && cp {input.r2} {output.r2}) > {log} 2>&1"
