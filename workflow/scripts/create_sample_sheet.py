import os
import re
import sys

## write to log file
sys.stderr = open(snakemake.log[0], "w")

inpath = snakemake.params.inpath
renaming = snakemake.params.renaming
sample_csv = snakemake.input[0]


def rename_fastqs(path):
    samples = []
    
    fastqs = [file for file in os.listdir(path) if file.endswith(".fastq.gz")]
    if not fastqs:
        print(
            f"Error: There are no fastq files in the directory. Have you used the correct path: {path}?"
        )
        raise Exception(
            f"There are no fastq files in the directory. Have you used the correct path: {path}?"
        )

    if renaming:
        print(
            "Renaming fastq files, e.g. from sampleID_S40_L001_R1_001.fastq.gz to sampleID_R1.fastq.gz"
        )
    else:
        print("Fastq files will not be renamed")

    for fastq in fastqs:
        ## renaming from e.g. sampleID_S40_L001_R1_001.fastq.gz to sampleID_R1.fastq.gz
        fastq_new = re.sub(r"_S\d{0,2}_L001", "", fastq)
        fastq_new = re.sub(r"_001.fastq", ".fastq", fastq_new)

        sample = (re.search("(.*)_R[1-2].fastq.gz", fastq_new)).group(1)
        if sample not in samples and sample != "Undetermined":
            samples.append(sample)

        if renaming:
            os.system(f"mv {path} {fastq} {path} {fastq_new}")

    return samples


def write_sample_sheet(samples, path, outfile):
    # os.system(f"touch {outfile}")

    with open(outfile, "w") as sheet:
        sheet.write("sample_name,fq1,fq2\n")

        for sample in samples:
            sheet.write(
                f"{sample},{path} {sample}_R1.fastq.gz,{path} {sample}_R2.fastq.gz\n"
            )


samples = rename_fastqs(inpath)
write_sample_sheet(samples, inpath, sample_csv)
