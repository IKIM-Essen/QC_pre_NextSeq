import pandas as pd
import altair as alt
import json
import re
import os
import sys

## write to log file
sys.stderr = open(snakemake.log[0], "w")

## input files
stat_files = snakemake.input.stats
bracken_domain = snakemake.input.bracken
json_files = snakemake.input.jsons


## output files
contamination_html = snakemake.output.human_cont_html  # "conta.html" #
domain_abundance_html = snakemake.output.domain_abd_html  #'domain_abundance.html' #
filtering_html = snakemake.output.read_summary_html  #'filtering_summary.html' #
summary_out_csv = snakemake.output.summary_csv  # "summary.csv" #

## variables
color_red = "#e03e3e"
color_green = "#6aa84f"


## functions
def get_human_contamination_df(stat_files):
    sum_dict = {}
    for stats_path in stat_files:
        # file=stats_path.split("/")[-1]
        sample = (re.search("(.*)_stats.txt", os.path.basename(stats_path))).group(1)
        sample_sum_dict = {}
        with open(stats_path, "r") as stats:
            for line in stats:
                if line.startswith("SN	sequences"):
                    total = line.split(":")[-1].strip()
                    continue
                elif line.startswith("SN	reads mapped"):
                    mapped = line.split(":")[-1].strip()
                    # sample_sum_dict["reads_mapped"] =int(mapped)

                    prc = int(mapped) / int(total)
                    sample_sum_dict["Human"] = prc  # "%.6f" %
                    break

        sum_dict[sample] = sample_sum_dict

    human_cont_df = pd.DataFrame.from_dict(sum_dict, orient="index")
    human_cont_df = human_cont_df.reset_index()
    human_cont_df.rename(columns={"index": "sample"}, inplace=True)
    human_cont_df.sort_values(by=["sample"], inplace=True)

    # human_cont_df.to_csv(out_csv)
    return human_cont_df


def plot_human_contamination(human_cont_df, out_html):
    slider = alt.binding_range(
        min=0, max=100, step=0.5, name="max human contamination:"
    )
    # selector = alt.param(name='SelectorName', value=50, bind=slider)
    selector = alt.selection_point(
        name="SelectorName", fields=["max_contamination"], bind=slider, value=50
    )

    # ,title="% human contamination"
    base_chart = (
        alt.Chart(human_cont_df)
        .encode(
            alt.X("Human:Q")
            .axis(format="%", labelFontSize=12, titleFontSize=15)
            .title("human contamination"),
            alt.Y("sample:N").axis(labelFontSize=12, titleFontSize=15),
        )
        .add_params(selector)
        .properties(width="container")
        .interactive()
    )

    bars = base_chart.mark_bar().encode(
        color=alt.condition(
            (alt.datum.Human * 100) >= selector.max_contamination,
            alt.value(color_red),
            alt.value(color_green),
        )
    )

    chart_text = base_chart.mark_text(
        align="center",
        baseline="middle",
        dx=20,
        fontSize=12,
    ).encode(
        text=alt.Text("Human:Q", format=".2%"),
    )

    full_chart = bars + chart_text
    full_chart.save(out_html)


def get_domain_abundance_df(infile):
    df = pd.read_table(infile)
    # only keep bracken fraction columns & species name
    df = df[df.columns.drop(list(df.filter(regex="bracken_num$|^taxonomy")))]

    df.columns = df.columns.str.replace(".bracken_frac", "", regex=False)

    # transpose df, change index & column names to get format for plotting
    df_trans = df.transpose()
    df_trans.reset_index(inplace=True)

    df_trans["index"][0] = "sample"
    df_trans.columns = df_trans.iloc[0]
    df_trans.drop(df_trans.index[0], inplace=True)

    df_trans.sort_values(by=["sample"], inplace=True)

    return df_trans


def plot_domain_abundance(domain_abundance_df, out_html):
    melt_df = domain_abundance_df.melt(
        id_vars=["sample"], var_name="Domain", value_name="share"
    )

    bars = (
        alt.Chart(melt_df, title="Relative abundance of domains")
        .mark_bar()
        .transform_calculate(
            combined_tooltip="datum.Domain + ': ' + format(datum.share, '.2%')"
        )
        .encode(
            alt.X("sample:N").axis(labelFontSize=12, titleFontSize=15).title("Sample"),
            alt.Y("sum(share)", stack="normalize")
            .axis(format="%", labelFontSize=12, titleFontSize=15)
            .title("Relative abundance"),
            color=alt.Color("Domain"),
            tooltip="combined_tooltip:N",
        )
        .properties(width="container", height=600)
    )

    bars = bars.configure_legend(
        titleFontSize=15, labelFontSize=12, labelFontStyle="italic"
    ).configure_title(fontSize=18)

    bars.save(out_html)


def get_qc_filtering_dataframes(json_files):
    ## contains number of reads before and after filtering & number of bases
    filtering_results_dict = {}
    ## contains number of reads after filtering & rate of Q30 bases
    read_quality_dict = {}

    for jsonfile in json_files:
        sample_filt_results_dict = {}
        sample_read_quality_dict = {}

        sample = (re.search("(.*).fastp.json", os.path.basename(jsonfile))).group(1)

        with open(jsonfile, "r") as f:
            fastp = json.load(f)

        total_reads = fastp["summary"]["after_filtering"]["total_reads"]

        # save number of bases as Mbp
        bases = fastp["summary"]["after_filtering"]["total_bases"]
        sample_filt_results_dict["total_bases"] = "{} Mbp".format(
            round((bases / 1000000))
        )

        sample_filt_results_dict["before filtering"] = fastp["summary"][
            "before_filtering"
        ]["total_reads"]
        sample_filt_results_dict["after filtering"] = total_reads

        filtering_results_dict[sample] = sample_filt_results_dict

        sample_read_quality_dict["Total reads"] = total_reads

        q_30 = fastp["summary"]["after_filtering"]["q30_rate"]
        sample_read_quality_dict["Q30 bp (%)"] = round((q_30 * 100), 3)

        read_quality_dict[sample] = sample_read_quality_dict

    filtering_results_df = pd.DataFrame.from_dict(
        filtering_results_dict, orient="index"
    )
    filtering_results_df = filtering_results_df.reset_index()
    filtering_results_df.rename(columns={"index": "sample"}, inplace=True)
    filtering_results_df.sort_values(
        by=["after filtering"], inplace=True, ignore_index=True
    )

    read_quality_df = pd.DataFrame.from_dict(read_quality_dict, orient="index")
    read_quality_df.index.name = "sample"
    read_quality_df.sort_index(inplace=True)

    return filtering_results_df, read_quality_df


def plot_filtering_results(filt_results_df, out_html):
    no_bases_df = filt_results_df[filt_results_df.columns.drop(["total_bases"])]
    melt_df = no_bases_df.melt(
        id_vars=["sample"], var_name="Status", value_name="number"
    )

    melt_df["total_bases"] = filt_results_df["total_bases"]
    melt_df = melt_df.fillna("")

    bars = (
        alt.Chart(melt_df)
        .mark_bar()
        .transform_calculate(
            combined_tooltip="datum.Status + ': ' + format(datum.number, ',')"
        )
        .encode(
            alt.Y("sample:N")
            .sort("-x")
            .axis(labelFontSize=12, titleFontSize=15)
            .title("Sample"),
            alt.X("number", stack=None)
            .axis(labelFontSize=12, titleFontSize=15)
            .title("Number of reads"),
            color=alt.Color("Status")
            .scale(range=[color_green, color_red])
            .legend(titleFontSize=15, labelFontSize=12),
            tooltip="combined_tooltip:N",
        )
        .properties(width="container")
    )

    chart_text = bars.mark_text(
        align="center",
        baseline="middle",
        dx=25,
        fontSize=12,
    ).encode(
        text=alt.Text("total_bases"),
        color=alt.value("black"),
    )

    full_chart = bars + chart_text
    full_chart.save(out_html)


def save_summary_csv(domain_abundance_df, human_cont_df, read_quality_df, outfile):
    domain_abundance_for_csv = domain_abundance_df.copy()
    domain_abundance_for_csv.set_index("sample", inplace=True)

    human_cont_for_csv = human_cont_df.copy()
    human_cont_for_csv.set_index("sample", inplace=True)

    df_all_for_csv = pd.concat([domain_abundance_for_csv, human_cont_for_csv], axis=1)

    # header = ["Human", "Bacteria", "Eukaryota", "Archaea", "Viruses"]
    # new_cols = [s + " (%)" for s in header]
    # df_all_for_csv = df_all_for_csv[header]
    # df_all_for_csv.columns = new_cols

    df_all_for_csv = df_all_for_csv.mul(100)
    df_all_for_csv = df_all_for_csv.astype("float64").round(3)

    df_all_for_csv = pd.concat([read_quality_df, df_all_for_csv], axis=1)

    df_all_for_csv.to_csv(outfile)


## running
human_cont_df = get_human_contamination_df(stat_files)
plot_human_contamination(human_cont_df, contamination_html)

domain_abundance_df = get_domain_abundance_df(bracken_domain)
plot_domain_abundance(domain_abundance_df, domain_abundance_html)

filtering_results_df, read_quality_df = get_qc_filtering_dataframes(json_files)
plot_filtering_results(filtering_results_df, filtering_html)

save_summary_csv(domain_abundance_df, human_cont_df, read_quality_df, summary_out_csv)
