import os
import re
path="/groups/ds/metagenomes/231218_Miseq/"
outfile="config/pep/samples_231218.csv"


def rename_fastqs(path):
    fastqs=os.listdir(path)
    samples=[]
    for fastq in fastqs:
        fastq_new = re.sub(r"_S\d{0,2}_L001", "", fastq)
        fastq_new = re.sub(r"_001.fastq", ".fastq", fastq_new)
        sample=fastq_new.split("_")[0]
        if sample not in samples:
            samples.append(sample)
        os.system(f"mv {path}{fastq} {path}{fastq_new}")
    return(samples)

def write_sample_sheet(samples,outfile):
    os.system(f"touch {outfile}")
    with open(outfile,"w") as sheet:
        sheet.write("sample_name,fq1,fq2\n")
        for sample in samples:
            sheet.write(f"{sample},{path}{sample}_R1.fastq.gz,{path}{sample}_R2.fastq.gz\n")

samples=rename_fastqs(path)
write_sample_sheet(samples,outfile)
