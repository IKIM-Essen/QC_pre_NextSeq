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

        match = re.search("(.*)_R[1-2].fastq.gz", fastq_new)
        if not match:
            continue

        sample = match.group(1)
        if sample not in samples and sample != "Undetermined":
            samples.append(sample)

        if renaming:
            os.system(f"mv {path}{fastq} {path}{fastq_new}")

    return samples


def write_sample_sheet(samples, path, outfile):
    with open(outfile, "w") as sheet:
        sheet.write("sample_name,fq1,fq2\n")

        for sample in samples:
            r1 = f"{path}{sample}_R1.fastq.gz"
            r2 = f"{path}{sample}_R2.fastq.gz"

            # check if FASTQs existieren and are not empty
            if not os.path.exists(r1) or not os.path.exists(r2):
                print(f"Skipping sample {sample}: FASTQ files missing.")
                continue

            if os.path.getsize(r1) == 0 or os.path.getsize(r2) == 0:
                print(f"Skipping sample {sample}: FASTQ files are empty.")
                continue

            sheet.write(f"{sample},{r1},{r2}\n")


samples = rename_fastqs(inpath)
write_sample_sheet(samples, inpath, sample_csv)
