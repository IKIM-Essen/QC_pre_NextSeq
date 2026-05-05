import os
import sys

sys.stderr = open(snakemake.log[0], "w")

valid_samples = []

for sig in snakemake.input.sigs:
    sample = os.path.basename(sig).replace(".sig", "")
    if os.path.exists(sig):
        valid_samples.append(sample)

os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)

with open(snakemake.output[0], "w") as f:
    for s in sorted(valid_samples):
        f.write(s + "\n")

print(f"[INFO] Wrote {len(valid_samples)} valid samples")


