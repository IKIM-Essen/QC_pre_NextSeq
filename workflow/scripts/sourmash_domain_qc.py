#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

# Snakemake übergibt automatisch die Pfade
gather_files = snakemake.input.gather
output_file = Path(snakemake.output[0])
log_file = Path(snakemake.log[0])

# Optional: Logfile zum Schreiben öffnen
with open(log_file, "w") as log:
    log.write(f"Processing {len(gather_files)} gather files...\n")

# Alle Gather-Dateien laden und in ein DataFrame zusammenführen
df_list = []
for f in gather_files:
    df = pd.read_csv(f, sep="\t")
    # Nur relevante Domains behalten: Bacteria, Viruses, Human
    df = df[df['name'].str.contains("gtdb|virus|human", case=False)]
    # Sample aus Dateiname extrahieren
    sample = Path(f).stem.replace("_gather", "")
    df['sample'] = sample
    df_list.append(df)

# Alle zusammenführen
merged_df = pd.concat(df_list, ignore_index=True)

# Optional: Pivot auf Domain-Level
domain_qc = merged_df.pivot_table(
    index='sample',
    columns='name',
    values='f_match',  # oder die Spalte, die den Anteil der Reads enthält
    fill_value=0
)

# Speichern
domain_qc.to_csv(output_file, sep="\t")

with open(log_file, "a") as log:
    log.write(f"Domain QC report written to {output_file}\n")

