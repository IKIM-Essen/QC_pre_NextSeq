import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


reports = snakemake.input.reports
outfile = snakemake.output.csv

tax_ids = {
    "bacteria": 2,
    "human": 9606,
    "archaea": 2157,
    "virus": 10239,
    "unclassified": 0,
}

header_ls = ["prct", "total_reads", "lvl_reads", "lvl", "tax_id", "name"]
results_dict = {}

for report in reports:
    report_df = pd.read_table(report, names=header_ls)
    sample_results_dict = {}
    for ref in tax_ids:
        line = report_df.loc[report_df["tax_id"] == tax_ids[ref]]
        if not line.empty:
            if line.iloc[0]["lvl"] == "S":
                no_reads = line.iloc[0]["lvl_reads"]
            else:
                no_reads = line.iloc[0]["total_reads"]
            prct_reads = line.iloc[0]["prct"]

        else:
            no_reads = 0
            prct_reads = 0

        if tax_ids[ref] == 0:
            no_uncl_reads = no_reads

        sample_results_dict[f"{ref} (%)"] = prct_reads
        sample_results_dict[f"{ref} (#)"] = f"{no_reads:,}"

    line = report_df.loc[report_df["tax_id"] == 1]
    total_reads = no_uncl_reads + line.iloc[0]["total_reads"]

    sample_results_dict["total reads (#)"] = f"{total_reads:,}"

    sample = report[report.rfind("/") + 1 : report.rfind("_report")]
    results_dict[sample] = sample_results_dict


results_df = pd.DataFrame.from_dict(results_dict, orient="index")
results_df.sort_index(inplace=True)
results_df.index.name = "sample"

# move total reads column to the front of the table
first_column = results_df.pop("total reads (#)")
results_df.insert(0, "total reads (#)", first_column)

results_df.to_csv(outfile)