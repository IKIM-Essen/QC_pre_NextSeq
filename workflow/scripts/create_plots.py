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
kaiju_inputs = list(snakemake.input.kaiju)
json_files = snakemake.input.jsons

## output files
contamination_html = snakemake.output.human_cont_html
domain_abundance_html = snakemake.output.domain_abd_html
filtering_html = snakemake.output.read_summary_html
summary_out_csv = snakemake.output.summary_csv
genus_abundance_html = snakemake.output.genus_abd_html
genus_top10_csv = snakemake.output.genus_top10_csv
# Machine-readable per-sample QC metrics for the sample registry (raw numbers,
# NOT display-formatted like filtering_summary.csv). Columns:
#   sample  raw_read_pairs  trimmed_read_pairs  q30_pct  human_pct
qc_metrics_tsv = snakemake.output.qc_metrics_tsv

## variables
color_red = "#e03e3e"
color_green = "#6aa84f"

DOMAIN_ORDER = ["Bacteria", "Eukaryota", "Archaea", "Viruses", "Unclassified"]
DOMAIN_COLORS = ["#4e79a7", "#f28e2b", "#e15759", "#59a14f", "#9d9d9d"]

GENUS_PALETTE = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac",
]


def _pick_input(inputs, patterns):
    for p in inputs:
        for pat in patterns:
            if pat in p:
                return p
    return inputs[0] if inputs else None


def get_human_contamination_df(stat_files):
    def _first_int(s):
        m = re.findall(r"\d+", s)
        return int(m[0]) if m else 0

    sum_dict = {}
    for stats_path in stat_files:
        sample = (re.search("(.*)_stats.txt", os.path.basename(stats_path))).group(1)
        sample_sum_dict = {}
        with open(stats_path, "r") as stats:
            total = None
            mapped = None
            for line in stats:
                if re.match(r"^SN\s+sequences:", line):
                    total = _first_int(line.split(":", 1)[-1])
                elif re.match(r"^SN\s+reads mapped:", line):
                    # anchor on the colon so this does NOT also match
                    # "SN  reads mapped and paired:" (a later line that would
                    # otherwise overwrite `mapped` with the wrong count)
                    mapped = _first_int(line.split(":", 1)[-1])

            if total and mapped is not None:
                prc = mapped / total
                sample_sum_dict["Human"] = prc
            else:
                sample_sum_dict["Human"] = 0.0

        sum_dict[sample] = sample_sum_dict

    human_cont_df = pd.DataFrame.from_dict(sum_dict, orient="index")
    human_cont_df = human_cont_df.reset_index()
    human_cont_df.rename(columns={"index": "sample"}, inplace=True)
    human_cont_df.sort_values(by=["sample"], inplace=True)
    return human_cont_df


def plot_human_contamination(human_cont_df, out_html):
    slider = alt.binding_range(
        min=0, max=100, step=0.5, name="max human contamination:"
    )
    selector = alt.selection_point(
        name="SelectorName", fields=["max_contamination"], bind=slider, value=50
    )

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

    (bars + chart_text).save(out_html)


def get_domain_abundance_df(infile):
    df = pd.read_table(infile)

    cols = {c.lower(): c for c in df.columns}
    file_col = cols.get("file")
    pct_col = cols.get("percent")
    name_col = cols.get("taxon_name", cols.get("name"))

    if file_col is None or pct_col is None or name_col is None:
        raise ValueError(f"Missing required columns. Found: {list(df.columns)}")

    # remove repeated header rows
    df = df[df[file_col] != "file"]

    df[pct_col] = pd.to_numeric(df[pct_col], errors="coerce")
    df = df[df[pct_col].notna()]

    df["sample"] = (
        df[file_col].astype(str)
        .str.replace(r"\.out$", "", regex=True)
        .str.replace(r".*/", "", regex=True)
    )
    df["share"] = df[pct_col] / 100.0

    # Domain from taxon_name like "Bacteria;Actinomycetota;"
    df["Domain"] = df[name_col].astype(str).str.split(";").str[0].str.strip()

    df["Domain"] = df["Domain"].replace(
        {
            "unclassified": "Unclassified",
            "cannot be assigned to a (non-viral) phylum": "Unclassified",
        }
    )


    out = (
        df.groupby(["sample", "Domain"], as_index=False)["share"]
        .sum()
        .pivot(index="sample", columns="Domain", values="share")
        .fillna(0)
        .reset_index()
    )
    out.columns.name = None
    return out


def plot_domain_abundance(domain_abundance_df, out_html):
    keep_cols = ["sample"] + [d for d in DOMAIN_ORDER if d in domain_abundance_df.columns]
    domain_abundance_df = domain_abundance_df[keep_cols]

    melt_df = domain_abundance_df.melt(
        id_vars=["sample"], var_name="Domain", value_name="share"
    )

    color_scale = alt.Scale(domain=DOMAIN_ORDER, range=DOMAIN_COLORS)

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
            color=alt.Color("Domain", scale=color_scale),
            tooltip="combined_tooltip:N",
        )
        .properties(width="container", height=600)
    )

    bars = bars.configure_legend(
        titleFontSize=15, labelFontSize=12, labelFontStyle="italic"
    ).configure_title(fontSize=18)

    bars.save(out_html)


def get_genus_top10_per_sample(infile):
    df = pd.read_table(infile)

    cols = {c.lower(): c for c in df.columns}
    file_col = cols.get("file")
    pct_col = cols.get("percent")
    name_col = cols.get("taxon_name", cols.get("name"))

    if file_col is None or pct_col is None or name_col is None:
        raise ValueError(f"Missing required columns. Found: {list(df.columns)}")

    df = df[df[file_col] != "file"]
    df[pct_col] = pd.to_numeric(df[pct_col], errors="coerce")
    df = df[df[pct_col].notna()]

    df["sample"] = (
        df[file_col].astype(str)
        .str.replace(r"\.out$", "", regex=True)
        .str.replace(r".*/", "", regex=True)
    )
    df["share"] = df[pct_col] / 100.0
    df["Genus"] = df[name_col].astype(str).str.strip()

    df["bucket"] = df["Genus"]
    df.loc[df["Genus"].str.contains("unclassified", case=False, na=False), "bucket"] = "Unclassified"
    df.loc[df["Genus"].str.contains("cannot be assigned", case=False, na=False), "bucket"] = "Unclassified"
    df = df[df["Genus"] != "Viruses"]

    classified = df[df["bucket"] != "Unclassified"].copy()
    classified["rank"] = classified.groupby("sample")["share"].rank(method="first", ascending=False)

    top10 = classified[classified["rank"] <= 10].copy()
    top10_table = (
        top10.sort_values(["sample", "share"], ascending=[True, False])
        .assign(percent=lambda x: (x["share"] * 100).round(2))
        [["sample", "Genus", "percent"]]
    )

    top10_keys = set(zip(top10["sample"], top10["Genus"]))

    def assign_bucket(row):
        if row["bucket"] == "Unclassified":
            return "Unclassified"
        if (row["sample"], row["Genus"]) in top10_keys:
            return row["Genus"]
        return "Rest"

    df["plot_genus"] = df.apply(assign_bucket, axis=1)
    plot_df = (
        df.groupby(["sample", "plot_genus"], as_index=False)["share"]
        .sum()
    )

    return plot_df, top10_table


def plot_genus_composition(plot_df, out_html):
    genus_order = sorted(
        [g for g in plot_df["plot_genus"].unique() if g not in {"Rest", "Unclassified"}]
    ) + ["Rest", "Unclassified"]

    color_range = GENUS_PALETTE[: max(0, len(genus_order) - 2)] + ["#c7c7c7", "#7f7f7f"]

    bars = (
        alt.Chart(plot_df, title="Genus composition (Top 10 per sample)")
        .mark_bar()
        .encode(
            alt.X("sample:N").axis(labelFontSize=12, titleFontSize=15).title("Sample"),
            alt.Y("sum(share):Q", stack="normalize")
            .axis(format="%", labelFontSize=12, titleFontSize=15)
            .title("Relative abundance"),
            color=alt.Color(
                "plot_genus:N",
                sort=genus_order,
                scale=alt.Scale(domain=genus_order, range=color_range),
            ),
            tooltip=[
                alt.Tooltip("sample:N"),
                alt.Tooltip("plot_genus:N", title="Genus"),
                alt.Tooltip("share:Q", format=".2%"),
            ],
        )
        .properties(width="container", height=600)
    )

    bars.save(out_html)


def write_top10_table(df_top10, out_csv):
    df_out = df_top10.copy()
    df_out = df_out.sort_values(["sample", "percent"], ascending=[True, False]).reset_index(drop=True)
    df_out["rank"] = df_out.groupby("sample").cumcount() + 1
    df_out = df_out[["sample", "rank", "Genus", "percent"]]
    df_out.columns = ["sample", "rank", "genus", "percent"]
    df_out["percent"] = df_out["percent"].round(3)
    df_out.to_csv(out_csv, index=False)


def get_qc_filtering_dataframes(json_files):
    filtering_results_dict = {}
    read_quality_dict = {}

    for jsonfile in json_files:
        sample_filt_results_dict = {}
        sample_read_quality_dict = {}

        sample = (re.search("(.*).fastp.json", os.path.basename(jsonfile))).group(1)

        with open(jsonfile, "r") as f:
            fastp = json.load(f)

        total_reads = fastp["summary"]["after_filtering"]["total_reads"]

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

    (bars + chart_text).save(out_html)


def save_summary_csv(domain_abundance_df, human_cont_df, read_quality_df, outfile):
    domain_abundance_for_csv = domain_abundance_df.copy()
    domain_abundance_for_csv.set_index("sample", inplace=True)

    human_cont_for_csv = human_cont_df.copy()
    human_cont_for_csv.set_index("sample", inplace=True)

    df_all_for_csv = pd.concat(
        [read_quality_df, domain_abundance_for_csv, human_cont_for_csv], axis=1
    )

    ordered_cols = [
        "Total reads",
        "Q30 bp (%)",
        "Human",
        "Bacteria",
        "Eukaryota",
        "Archaea",
        "Viruses",
        "Unclassified",
    ]
    df_all_for_csv = df_all_for_csv.reindex(columns=ordered_cols, fill_value=0)

    df_all_for_csv["Total reads"] = (
        pd.to_numeric(df_all_for_csv["Total reads"], errors="coerce")
        .fillna(0)
        .map(lambda x: f"{int(x):,}")
    )

    df_all_for_csv["Q30 bp (%)"] = (
        pd.to_numeric(df_all_for_csv["Q30 bp (%)"], errors="coerce")
        .fillna(0)
        .map(lambda x: f"{x:.2f}")
    )

    percent_cols = ["Human", "Bacteria", "Eukaryota", "Archaea", "Viruses", "Unclassified"]
    for col in percent_cols:
        df_all_for_csv[col] = (
            pd.to_numeric(df_all_for_csv[col], errors="coerce")
            .fillna(0)
            .mul(100)
            .map(lambda x: f"{x:.2f}")
        )

    df_all_for_csv.rename(
        columns={
            "Human": "Human (%)",
            "Bacteria": "Bacteria (%)",
            "Eukaryota": "Eukaryota (%)",
            "Archaea": "Archaea (%)",
            "Viruses": "Viruses (%)",
            "Unclassified": "Unclassified (%)",
        },
        inplace=True,
    )

    df_all_for_csv.to_csv(outfile)


## running
human_cont_df = get_human_contamination_df(stat_files)
plot_human_contamination(human_cont_df, contamination_html)

# kaiju is optional. In QC mode there are no kaiju inputs, but the Snakemake
# rules still declare the domain/genus/top10 outputs, so we MUST create those
# files (empty placeholders) to satisfy the DAG. domain_abundance_df stays empty
# so the summary CSV simply has no domain columns.
if kaiju_inputs:
    kaiju_domain_file = _pick_input(kaiju_inputs, ["merged.kaiju_domain", "merged.kaiju_phylum"])
    domain_abundance_df = get_domain_abundance_df(kaiju_domain_file)
    plot_domain_abundance(domain_abundance_df, domain_abundance_html)
    if genus_abundance_html:
        kaiju_genus_file = _pick_input(kaiju_inputs, ["merged.kaiju_genus"])
        if kaiju_genus_file:
            genus_plot_df, genus_top10_df = get_genus_top10_per_sample(kaiju_genus_file)
            plot_genus_composition(genus_plot_df, genus_abundance_html)
            write_top10_table(genus_top10_df, genus_top10_csv)
        else:
            raise ValueError("genus_abundance_html requested but no merged.kaiju_genus input found")
else:
    # QC mode: write empty placeholders so all declared outputs exist.
    domain_abundance_df = pd.DataFrame(columns=["sample"])
    _placeholder = ("<html><body><p>kaiju taxonomic profiling was not run in "
                    "QC mode. Run the pipeline with mode=diversity for domain/"
                    "genus abundance.</p></body></html>")
    if domain_abundance_html:
        with open(domain_abundance_html, "w") as fh:
            fh.write(_placeholder)
    if genus_abundance_html:
        with open(genus_abundance_html, "w") as fh:
            fh.write(_placeholder)
    if genus_top10_csv:
        # minimal valid CSV so genus_top10_report (rbt csv-report) doesn't choke
        with open(genus_top10_csv, "w") as fh:
            fh.write("sample,rank,genus,percent\n")


filtering_results_df, read_quality_df = get_qc_filtering_dataframes(json_files)
plot_filtering_results(filtering_results_df, filtering_html)

save_summary_csv(domain_abundance_df, human_cont_df, read_quality_df, summary_out_csv)


def save_qc_metrics_tsv(filt_df, read_quality_df, human_cont_df, outfile):
    """Machine-readable per-sample QC metrics for the sample registry: RAW
    numbers (no thousands separators, no % strings), one row per sample.
    fastp counts reads (R1+R2); the registry tracks PAIRS, so before/after
    filtering counts are halved. Human contamination is mapped/total (a
    fraction) turned into a percent. Missing values are left blank."""
    import csv as _csv

    filt = filt_df.set_index("sample")
    human = human_cont_df.set_index("sample")
    rq = read_quality_df  # index is already 'sample'

    def _pairs(reads):
        try:
            return int(round(float(reads) / 2.0))
        except (TypeError, ValueError):
            return ""

    def _num(v, ndigits=None):
        try:
            f = float(v)
            return round(f, ndigits) if ndigits is not None else f
        except (TypeError, ValueError):
            return ""

    samples = sorted(set(filt.index) | set(rq.index) | set(human.index))
    with open(outfile, "w", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "raw_read_pairs", "trimmed_read_pairs",
                    "q30_pct", "human_pct"])
        for s in samples:
            raw = _pairs(filt.at[s, "before filtering"]) if s in filt.index else ""
            trimmed = _pairs(filt.at[s, "after filtering"]) if s in filt.index else ""
            q30 = _num(rq.at[s, "Q30 bp (%)"], 3) if s in rq.index else ""
            hum = _num(human.at[s, "Human"]) if s in human.index else ""
            hum = round(hum * 100, 4) if hum != "" else ""
            w.writerow([s, raw, trimmed, q30, hum])


save_qc_metrics_tsv(filtering_results_df, read_quality_df, human_cont_df, qc_metrics_tsv)
