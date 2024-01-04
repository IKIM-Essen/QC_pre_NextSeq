import numpy as np
import pandas as pd
import altair as alt
from distinctipy import distinctipy

sys.stderr = open(snakemake.log[0], "w")

infile = snakemake.input[0]
threshold = snakemake.params.threshold
level = snakemake.wildcards.level
out_html = snakemake.output[0]


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

    # replace all NaN with 0, because NaN can cause plotting problems
    df_plot.replace(np.nan, 0, inplace=True)
    # sort samples by name
    df_plot.sort_values(by=["sample"], inplace=True)

    return df_plot


def get_plot_title(level):
    level_to_title = {
        "genus": "genera",
        "family": "families",
        "class": "classes",
        "phylum": "phyla",
    }

    level_title = "Relative abundance of bacterial {}".format(level_to_title[level])
    return level_title


def get_hex_color_list(bracken_df):
    # get a list of distinct colors
    colors = distinctipy.get_colors(len(bracken_df.columns))
    hex_clrs=[]
    for color in colors:
        hex_clrs.append(distinctipy.get_hex(color))
    return(hex_clrs)


# stacked bar plot with relative abundance of species reads per sample
def plot_bracken(bracken_df, level, out_html):
    # creates a dataframe column (level value) gets one row per sample
    melt_df = bracken_df.melt(id_vars=["sample"], var_name=level, value_name="share")

    # get a list of distinct colors in hex code
    clrs = get_hex_color_list(bracken_df)
    
    # get title for the plot level dependent
    level_title = get_plot_title(level)

    # tooltips displayed when hovering over the bars
    tooltip = f"datum.{level} + ': ' + format(datum.share, '.2%')"

    # stacked bar plot
    bars = (
        alt.Chart(melt_df, title=level_title)
        .mark_bar()
        .transform_calculate(combined_tooltip=tooltip)
        .encode(
            alt.X("sample:N").axis(labelFontSize=12, titleFontSize=15).title("Sample"),
            alt.Y("sum(share)", stack="normalize")
            .axis(format="%", labelFontSize=12, titleFontSize=15)
            .title("Relative abundance"),
            color=alt.Color(level).scale(range=clrs),
            tooltip="combined_tooltip:N",
        )
    )
    # configuring legend, bacterial names are written italic
    bars = bars.configure_legend(
        titleFontSize=15, labelFontSize=12, labelFontStyle="italic"
    ).configure_title(fontSize=18)

    # saving the plot as html file
    bars.save(out_html)


bracken_df = bracken2df(infile, threshold)
plot_bracken(bracken_df, level, out_html)
