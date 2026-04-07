import os
import re
import sys
from collections import defaultdict

# Write to log file
sys.stderr = open(snakemake.log[0], "w")

inpath = os.path.abspath(snakemake.params.inpath)
renaming = snakemake.params.renaming
sample_csv = snakemake.output[0]


FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)"
    r"(?:_S\d+_L\d{3})?"
    r"_R(?P<read>[12])"
    r"(?:_001)?"
    r"\.(?:fastq|fq)\.gz$"
)


def discover_fastqs(path):
    if not os.path.isdir(path):
        raise Exception(f"Input path is not a directory: {path}")

    files = sorted(
        [f for f in os.listdir(path) if f.endswith(".fastq.gz") or f.endswith(".fq.gz")]
    )
    if not files:
        raise Exception(f"No FASTQ files found in: {path}")

    pairs = defaultdict(dict) 

    for fname in files:
        m = FASTQ_RE.match(fname)
        if not m:
            continue

        sample = m.group("sample")
        read = m.group("read")

        if sample == "Undetermined":
            continue

        full = os.path.join(path, fname)

        # prefer name if both exist (e.g. _R1.fastq.gz over _R1_001.fastq.gz)
        current = pairs[sample].get(read)
        if current is None:
            pairs[sample][read] = full
        else:
            canon_current = bool(re.search(r"_R[12]\.(fastq|fq)\.gz$", os.path.basename(current)))
            canon_new = bool(re.search(r"_R[12]\.(fastq|fq)\.gz$", fname))
            if canon_new and not canon_current:
                pairs[sample][read] = full

    return pairs, files


def maybe_rename_files(path, files):
    if not renaming:
        print("Fastq files will not be renamed")
        return

    print("Renaming FASTQ files to canonical format sample_R1/2.fastq.gz where possible")
    for old in files:
        m = FASTQ_RE.match(old)
        if not m:
            continue

        sample = m.group("sample")
        read = m.group("read")
        new = f"{sample}_R{read}.fastq.gz"

        if old == new:
            continue

        src = os.path.join(path, old)
        dst = os.path.join(path, new)

        # avoid overwrite
        if os.path.exists(dst):
            print(f"Skip rename (target exists): {old} -> {new}")
            continue

        os.rename(src, dst)


def write_sample_sheet(pairs, outfile):
    with open(outfile, "w") as sheet:
        sheet.write("sample_name,fq1,fq2\n")

        for sample in sorted(pairs.keys()):
            r1 = pairs[sample].get("1")
            r2 = pairs[sample].get("2")

            if not r1 or not r2:
                print(f"Skipping sample {sample}: FASTQ pair incomplete (R1 or R2 missing).")
                continue

            if os.path.getsize(r1) == 0 or os.path.getsize(r2) == 0:
                print(f"Skipping sample {sample}: FASTQ files are empty.")
                continue

            sheet.write(f"{sample},{r1},{r2}\n")

maybe_rename_files(inpath, os.listdir(inpath))
pairs, _ = discover_fastqs(inpath)
write_sample_sheet(pairs, sample_csv)
