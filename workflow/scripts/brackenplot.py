import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from distinctipy import distinctipy

sys.stderr = open(snakemake.log[0], "w")

infile = snakemake.input[0]
threshold = snakemake.params.threshold
level=snakemake.wildcards.level
png = snakemake.output[0]


# reads in bracken output
def bracken2df(infile, threshold):
    df = pd.read_table(infile)
    # only keep bracken fraction columns & species name
    df = df[df.columns.drop(list(df.filter(regex="bracken_num$|^taxonomy")))]

    # all fractions below threshold summed up to a value for 'other' per sample
    sample_ls = list(df.filter(regex="bracken_frac"))
    other_row = ["other"]
    for sample in sample_ls:
        other_row.append(df.loc[df[sample] < threshold, sample].sum())
        # values below threshold set to NaN
        df.loc[df[sample] < threshold, sample] = np.nan
    df.sort_values(["name"], inplace=True)
    df.loc[-1] = other_row

    # rows with NaN for all samples are removed
    df = df.dropna(subset=sample_ls, how="all")

    df.columns = df.columns.str.replace(".bracken_frac", "", regex=False)

    column_ls = list(df.columns)
    column_ls[0] = "sample"

    # transpose df, change index & column names to get format for plotting
    df_plot = df.transpose()
    df_plot.reset_index(inplace=True)
    df_plot["index"][0] = "sample"
    df_plot.columns = df_plot.iloc[0]
    df_plot.drop(df_plot.index[0], inplace=True)

    return df_plot


# stacked bar plot with relative abundance of species reads per sample
def plot_bracken(bracken_df, png_file):
    plt.style.use("default")
    labelsize = 16
    colors = distinctipy.get_colors(len(bracken_df.columns))
    clrs = distinctipy.get_colormap(colors)

    ax = bracken_df.plot(
        x="sample",
        kind="bar",
        stacked=True,
        figsize=(12, 8),
        fontsize=(labelsize - 2),
        colormap=clrs,  # "gist_rainbow",
    )

    # ax.set_title('new Stacked Bar Graph',fontsize=20)

    # format axis labels & ticks
    ax.xaxis.label.set_size(labelsize)
    ax.set_xticklabels(bracken_df["sample"], rotation=45)
    ax.set_ylabel("relative abundance of species reads", fontsize=labelsize)
    # ax.yaxis.set_major_formatter(mtick.PercentFormatter())

    # remove right frame & set y-axis limit to 1
    ax.spines["right"].set_visible(False)
    # ax.spines['top'].set_visible(False)
    plt.ylim([0, 1])

    # format legend
    handles, labels = ax.get_legend_handles_labels()
    # species names in italic
    new_labels = []
    for label in labels:
        if label == "other":
            new_labels.append(label)
        else:
            new_labels.append("$\\it{" + label.replace(" ", "\ ") + "}$")

    legend = plt.legend(
        handles=handles,
        labels=new_labels,
        title=level,
        bbox_to_anchor=(1.35, 0.5),
        loc="center right",
        fontsize=(labelsize - 4),
    )
    legend.get_title().set_fontsize(labelsize - 2)

    # annotate percentage values (over 10%) to bars
    for p in ax.patches:
        if p.get_height() >= 0.1:
            perc = str(np.round((p.get_height() * 100), 2)) + "%"
            xpos = p.get_x() + (p.get_width() * 0.5)
            ypos = p.get_y() + (p.get_height() * 0.5)
            ax.annotate(
                perc,
                (xpos, ypos),
                fontsize=(labelsize - 8),
                ha="center",
                va="center",
                # rotation=45,
            )

    # save plot as png to output file
    plt.savefig(png_file, dpi=300, bbox_inches="tight")


df_plot = bracken2df(infile, threshold)
plot_bracken(df_plot, png)