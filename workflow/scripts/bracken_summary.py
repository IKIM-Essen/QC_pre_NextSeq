import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")


infile = snakemake.input.merged
outfile = snakemake.output.csv


def bracken2df(infile):
    df = pd.read_table(infile)
    # only keep bracken fraction columns & species name
    df = df[df.columns.drop(list(df.filter(regex="bracken_num$|^taxonomy")))]

    df.columns = df.columns.str.replace(".bracken_frac", "", regex=False)

    # transpose df, change index & column names to get format for plotting
    bracken_df = df.transpose()
    bracken_df.reset_index(inplace=True)

    new_cols=[s + " (%)" for s in bracken_df.iloc[0]]
    new_cols[0] = "sample"
    bracken_df.columns = new_cols

    bracken_df.set_index("sample", inplace=True)
    bracken_df.drop(bracken_df.index[0], inplace=True)

    bracken_df = bracken_df.mul(100)
    bracken_df = bracken_df.astype("float64").round(2)

    # sort samples by name
    bracken_df.sort_values(by=["sample"], inplace=True)

    return bracken_df

bracken_df = bracken2df(infile)
bracken_df.to_csv(outfile)
