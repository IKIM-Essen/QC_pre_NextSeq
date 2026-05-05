#!/usr/bin/env python3

import os
import re
import argparse
import yaml
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create sample sheet from FASTQ files (Illumina + ONT)"
    )
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Path to config.yaml (default: config/config.yaml)"
    )
    parser.add_argument(
        "-o", "--output",
        default="config/pep/samples.csv",
        help="Output sample sheet"
    )
    return parser.parse_args()


def load_config(config_path):
    with open(config_path) as f:
        return yaml.safe_load(f)


def find_fastqs(path):
    # Alle fastq.gz außer Undetermined_* Dateien
    return [f for f in os.listdir(path) if f.endswith(".fastq.gz") and not f.startswith("Undetermined")]

def normalize_illumina_name(fname):
    fname = re.sub(r"_S\d+_L\d{3}", "", fname)
    fname = re.sub(r"_001\.fastq\.gz$", ".fastq.gz", fname)
    return fname

def collect_samples(fastqs, path):
    samples = defaultdict(lambda: {"R1": None, "R2": None})

    for fq in fastqs:
        fq_path = os.path.join(path, fq)

        # Illumina paired-end: match bis _S#_R1/2_001.fastq.gz
        m = re.match(r"(.+)_S\d+_R([12])_001\.fastq\.gz", fq)
        if m:
            sample, read = m.groups()
            read_key = f"R{read}"
            # Nur für CSV normalisieren
            fq_path_csv = os.path.join(path, f"{sample}_R{read}.fastq.gz")
            samples[sample][read_key] = fq_path_csv
        else:
            # ONT oder Single-End Files
            sample = fq.replace(".fastq.gz", "")
            samples[sample]["R1"] = fq_path

    return samples


def write_sample_sheet(samples, outfile):
    with open(outfile, "w") as out:
        out.write("sample_name,fq1,fq2,technology\n")
        for sample in sorted(samples):
            fq1 = samples[sample]["R1"] or ""
            fq2 = samples[sample]["R2"] or ""
            tech = "illumina" if fq2 else "ont"
            out.write(f"{sample},{fq1},{fq2},{tech}\n")


def main():
    args = parse_args()
    
    # Pfad aus Config auslesen
    cfg = load_config(args.config)
    input_dir = cfg["sample-sheet"]["data-path"]

    fastqs = find_fastqs(input_dir)
    samples = collect_samples(fastqs, input_dir)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    write_sample_sheet(samples, args.output)

    print(f"Sample sheet written to {args.output}")


if __name__ == "__main__":
    main()


