import pandas as pd
import altair as alt
import json
import re
import os
import sys

# Redirect stderr to Snakemake log
sys.stderr = open(snakemake.log[0], "w")

# -------------------------------
# Inputs & Outputs
# -------------------------------
json_files = snakemake.input.jsons
stat_files = snakemake.input.stats
sourmash_tsv = snakemake.input.sourmash

human_cont_html = snakemake.output.human_cont_html
domain_abundance_html = snakemake.output.domain_abd_html
domain_blocks_html = snakemake.output.domain_blocks_html
read_summary_html = snakemake.output.read_summary_html
summary_out_csv = snakemake.output.summary_csv

color_red = "#e03e3e"
color_green = "#6aa84f"

# -------------------------------
# 1. Human Contamination Plot
# -------------------------------
def get_human_contamination_df(stat_files):
    sum_dict = {}
    for stats_path in stat_files:
        sample = re.search("(.*)_stats.txt", os.path.basename(stats_path)).group(1)
        sample_sum_dict = {}
        with open(stats_path) as f:
            for line in f:
                if line.startswith("SN\tsequences"):
                    total = int(line.split(":")[-1].strip())
                elif line.startswith("SN\treads mapped"):
                    mapped = int(line.split(":")[-1].strip())
                    sample_sum_dict["Human"] = mapped / total
                    break
        sum_dict[sample] = sample_sum_dict
    df = pd.DataFrame.from_dict(sum_dict, orient="index").reset_index()
    df.rename(columns={"index": "sample"}, inplace=True)
    df.sort_values("sample", inplace=True)
    return df

def plot_human_contamination(df, out_html):
    slider = alt.binding_range(min=0, max=100, step=0.5, name="max human contamination:")
    selector = alt.selection_point(name="SelectorName", fields=["max_contamination"], bind=slider, value=50)
    base_chart = alt.Chart(df).encode(
        alt.X("Human:Q", axis=alt.Axis(format="%", labelFontSize=12, titleFontSize=15), title="human contamination"),
        alt.Y("sample:N", axis=alt.Axis(labelFontSize=12, titleFontSize=15))
    ).add_params(selector).properties(width="container").interactive()

    bars = base_chart.mark_bar().encode(
        color=alt.condition((alt.datum.Human*100) >= selector.max_contamination,
                            alt.value(color_red), alt.value(color_green))
    )

    chart_text = base_chart.mark_text(align="center", baseline="middle", dx=20, fontSize=12).encode(
        text=alt.Text("Human:Q", format=".2%")
    )
    chart = bars + chart_text
    chart.save(out_html)

# -------------------------------
# 2. Domain Abundance & Blocks
# -------------------------------
def get_domain_abundance_df(sourmash_tsv):
    
    import pandas as pd
    import sys

    # Namedlist sicher abfangen
    if isinstance(sourmash_tsv, (list, tuple)):
        df = pd.concat(
            [pd.read_csv(f, sep="\t", low_memory=False) for f in sourmash_tsv],
            ignore_index=True
        )
    else:
        df = pd.read_csv(sourmash_tsv, sep="\t", low_memory=False)

    print("Sourmash columns:", df.columns.tolist(), file=sys.stderr)

    # 🔥 Sicherheitscheck (extrem empfohlen)
    REQUIRED_COLUMNS = ["filename", "name", "f_match", "query_name"]

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Sourmash TSV missing columns: {missing}")

    # ---------- Domain Mapping ----------
    def get_domain(row):
        src = (
            str(row.get("filename", "")) + " " +
            str(row.get("name", ""))
        ).lower()

        if "virus" in src:
            return "Viruses"
        elif "gtdb" in src or "bacter" in src:
            return "Bacteria"
        elif "human" in src or "homo sapiens" in src:
            return "Human"
        else:
            return "Other"

    # 🔥 DIESE ZEILE HAT DIR GEFEHLT
    df["Domain"] = df.apply(get_domain, axis=1)
    df["Domain"] = df["Domain"].astype(str)

    # ---------- Aggregation ----------
    domain = (
        df.groupby(["query_name", "Domain"])["f_match"]
        .sum()
        .unstack(fill_value=0)
    )

    # echte Prozentwerte erzeugen
    domain = domain.div(domain.sum(axis=1).replace(0, 1), axis=0)

    domain.reset_index(inplace=True)
    domain.rename(columns={"query_name": "sample"}, inplace=True)

    return domain


def plot_domain_abundance(df, out_html):
    
    melt_df = df.melt(
        id_vars="sample",
        var_name="Domain",
        value_name="share"
    )

    # Domains nach globaler Häufigkeit sortieren
    domain_order = (
        melt_df.groupby("Domain")["share"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )

    # Samples nach Bacteria Anteil sortieren
    if "Bacteria" in df.columns:
        sample_order = (
            df.sort_values("Bacteria", ascending=False)["sample"]
            .tolist()
        )
    else:
        sample_order = df["sample"].astype(str).tolist()
        melt_df["sample"] = melt_df["sample"].astype(str)

    color_scale = alt.Scale(
        domain=domain_order,
        range=[
            "#4CAF50",  # Bacteria - grün
            "#2196F3",  # Viruses - blau
            "#FFC107",  # Human - gelb
            "#9E9E9E"   # Other - grau
        ]
    )

    bars = (
        alt.Chart(melt_df,
                  title="Relative Domain Abundance per Sample")
        .mark_bar(size=28)
        .encode(
            x=alt.X(
                "sample:N",
                sort=sample_order,
                title="Sample",
                axis=alt.Axis(labelAngle=-30)
            ),

            y=alt.Y(
                "share:Q",
                title="Relative abundance",
                axis=alt.Axis(format="%")
            ),

            color=alt.Color(
                "Domain:N",
                scale=color_scale,
                legend=alt.Legend(title="Domain")
            ),

            tooltip=[
                alt.Tooltip("sample:N"),
                alt.Tooltip("Domain:N"),
                alt.Tooltip("share:Q", format=".2%")
            ]
        )
        .properties(
            width="container",
            height=520
        )
    )

    text = bars.mark_text(
        dy=-5,
        color="black"
    ).encode(
        text=alt.Text("share:Q", format=".0%")
    ).transform_filter(
        alt.datum.share > 0.05   # nur Labels >5%
    )

    domain_long = df.melt(
        id_vars="sample",
        var_name="Domain",
        value_name="Abundance"
    )

    chart = (
    alt.Chart(domain_long)
    .mark_bar()
    .encode(
        x="sample:N",
        y=alt.Y("Abundance:Q", stack="normalize"),
        color="Domain:N"
    )
)

    chart.save(out_html)

    print(domain_abundance_df.dtypes)
    print(domain_abundance_df["sample"].head())    

def plot_domain_blocks(df, out_html):
    melt_df = df.melt(id_vars=["sample"], var_name="Domain", value_name="share")
    melt_df = melt_df[melt_df["share"] > 0]  # Filter zero domains
    melt_df["Group"] = melt_df["Domain"].apply(lambda d: "Major domains" if d=="Bacteria" else "Minor domains")
    color_scale = alt.Scale(scheme="tableau20")

    major_chart = alt.Chart(melt_df[melt_df["Group"]=="Major domains"]).mark_bar(size=10, stroke="white", strokeWidth=0.5).encode(
        y=alt.Y("sample:N", title="Sample"),
        x=alt.X("sum(share):Q", stack="normalize", title="Relative abundance (%)"),
        color=alt.Color("Domain:N", scale=color_scale, title=None),
        tooltip=[alt.Tooltip("Domain:N"), alt.Tooltip("share:Q", format=".2%")]
    ).properties(width=800, height=80, title="Major domains")

    minor_chart = alt.Chart(melt_df[melt_df["Group"]=="Minor domains"]).mark_bar(size=10, stroke="white", strokeWidth=0.5).encode(
        y=alt.Y("sample:N", title="Sample"),
        x=alt.X("sum(share):Q", stack=None, title="Relative abundance (%)", scale=alt.Scale(domain=[0,1])),
        color=alt.Color("Domain:N", scale=color_scale, title=None),
        tooltip=[alt.Tooltip("Domain:N"), alt.Tooltip("share:Q", format=".2%")]
    ).properties(width=800, height=80, title="Minor domains (zoomed 0 – 1%)")

    chart = alt.vconcat(major_chart, minor_chart).resolve_scale(color="shared")
    chart.save(out_html)

# -------------------------------
# 3. QC Filtering Plots
# -------------------------------
def get_qc_filtering_dataframes(json_files):
    filtering_results_dict = {}
    read_quality_dict = {}
    for jsonfile in json_files:
        sample = re.search("(.*).fastp.json", os.path.basename(jsonfile)).group(1)
        with open(jsonfile) as f:
            fastp = json.load(f)
        total_reads = fastp["summary"]["after_filtering"]["total_reads"]
        bases = fastp["summary"]["after_filtering"]["total_bases"]
        filtering_results_dict[sample] = {
            "total_bases": f"{round(bases/1e6)} Mbp",
            "before filtering": fastp["summary"]["before_filtering"]["total_reads"],
            "after filtering": total_reads
        }
        read_quality_dict[sample] = {
            "Total reads": total_reads,
            "Q30 bp (%)": round(fastp["summary"]["after_filtering"]["q30_rate"]*100,3)
        }
    filtering_results_df = pd.DataFrame.from_dict(filtering_results_dict, orient="index").reset_index()
    filtering_results_df.rename(columns={"index":"sample"}, inplace=True)
    filtering_results_df.sort_values("after filtering", inplace=True, ignore_index=True)

    read_quality_df = pd.DataFrame.from_dict(read_quality_dict, orient="index")
    read_quality_df.index.name = "sample"
    read_quality_df.sort_index(inplace=True)

    return filtering_results_df, read_quality_df

def plot_filtering_results(filt_results_df, out_html):
    no_bases_df = filt_results_df[filt_results_df.columns.drop(["total_bases"])]
    melt_df = no_bases_df.melt(id_vars=["sample"], var_name="Status", value_name="number")
    melt_df["total_bases"] = filt_results_df["total_bases"]
    melt_df = melt_df.fillna("")
    bars = alt.Chart(melt_df).mark_bar().encode(
        alt.Y("sample:N").sort("-x").axis(labelFontSize=12, titleFontSize=15),
        alt.X("number", stack=None).axis(labelFontSize=12, titleFontSize=15),
        color=alt.Color("Status").scale(range=[color_green,color_red])
    ).properties(width="container")
    bars.save(out_html)

# -------------------------------
# Run everything
# -------------------------------
human_cont_df = get_human_contamination_df(stat_files)
plot_human_contamination(human_cont_df, human_cont_html)

domain_abundance_df = get_domain_abundance_df(snakemake.input.sourmash)
plot_domain_abundance(domain_abundance_df, domain_abundance_html)
plot_domain_blocks(domain_abundance_df, domain_blocks_html)

filtering_results_df, read_quality_df = get_qc_filtering_dataframes(json_files)
plot_filtering_results(filtering_results_df, read_summary_html)

# Save merged summary 
merged_df = pd.concat([domain_abundance_df.set_index("sample"), human_cont_df.set_index("sample")], axis=1)
merged_df.reset_index(inplace=True)
merged_df.to_csv(summary_out_csv, index=False)

