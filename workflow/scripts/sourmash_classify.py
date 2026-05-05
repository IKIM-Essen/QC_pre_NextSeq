import os
import subprocess

sig = snakemake.input.sig
out_csv = snakemake.output.summary
flag = snakemake.output.flag
log_file = snakemake.log[0]

os.makedirs(os.path.dirname(out_csv), exist_ok=True)
os.makedirs(os.path.dirname(flag), exist_ok=True)

cmd = [
    "sourmash", "classify",
    "--csv", out_csv,
    sig
]

with open(log_file, "w") as log:
    log.write("Running: " + " ".join(cmd) + "\n")
    subprocess.run(cmd, stdout=log, stderr=log, check=True)

# Flag schreiben, wenn erfolgreich
with open(flag, "w") as f:
    f.write("OK\n")


