import pandas as pd
import altair as alt
import re
import sys

sys.stderr = open(snakemake.log[0], "w")

stat_files=snakemake.input.stats
outfile=snakemake.output[0]

sum_dict={}
for stats_path in stat_files:
    file=stats_path.split("/")[-1]
    sample=(re.search('(.*)_stats.txt', file)).group(1)
    sample_sum_dict = {}
    with open(stats_path, "r") as stats:
        for line in stats:
            if line.startswith("SN	sequences"):
                total=line.split(":")[-1].strip()
                continue
            elif line.startswith("SN	reads mapped"):
                mapped=line.split(":")[-1].strip()
                sample_sum_dict["reads_mapped"] =int(mapped)

                prc=(int(mapped)/int(total))
                sample_sum_dict["contamination"] = "%.6f" % prc
                break

    sum_dict[sample]=sample_sum_dict

sum_df = pd.DataFrame.from_dict(sum_dict, orient="index")
sum_df = sum_df.reset_index()
sum_df.rename(columns={"index":"sample"}, inplace=True)


slider = alt.binding_range(min=0, max=100, step=0.5, name='max human contamination:')
#selector = alt.param(name='SelectorName', value=50, bind=slider)
selector = alt.selection_point(name="SelectorName", fields=['max_contamination'],
                                bind=slider, value=50)


base_chart=alt.Chart(sum_df,title="% human contamination").encode(
    alt.X('contamination:Q').axis(format='%',labelFontSize=12, titleFontSize=15).title('human contamination'),
    alt.Y('sample:N').axis(labelFontSize=12, titleFontSize=15),
).add_params(
   selector
).properties(width="container").interactive()

bars=base_chart.mark_bar().encode(color=alt.condition(
       (alt.datum.contamination*100) >= selector.max_contamination,
       alt.value('#9A0430'),
       alt.value('#6ea165')
   ))


chart_text = base_chart.mark_text(
    align='center',
    baseline='middle',
    dx=20,
    fontSize=12,
).encode(
    text=alt.Text("contamination:Q",format='.2%'),)


full_chart = bars + chart_text
full_chart.save(outfile)
