import os
import numpy
import pandas
import argparse
import collections
from tqdm import tqdm
from Bio import SeqIO
from pandas import ExcelWriter

parser = argparse.ArgumentParser()
parser.add_argument("--metadata", help="Enter Metadata")
parser.add_argument("--output", help="Enter output file in xlsx")
args = parser.parse_args()
metadata_url = args.metadata
output_url = args.output

metadata = pandas.read_csv(
    metadata_url, delimiter="\t", encoding="utf-8", low_memory=False
)
if len(metadata) > 0:
    index_keys = metadata["lineage"].unique().tolist()
    metadata["date"] = pandas.to_datetime(metadata["date"], format="%Y-%m-%d")
    month_values = sorted(metadata["date"])
    month_keys = [i.strftime("%b-%Y") for i in month_values]
    month_keys = pandas.DataFrame(month_keys)[0].unique().tolist()

    metadata = metadata.assign(collection_month=metadata["date"].dt.strftime("%b-%Y"))
    writer = ExcelWriter(output_url)

    all_voc_df_count = pandas.DataFrame(index=index_keys, columns=month_keys)
    all_voc_df_percent = pandas.DataFrame(index=index_keys, columns=month_keys)
    all_voc_df_combined = pandas.DataFrame(index=index_keys, columns=month_keys)

    for i in tqdm(month_keys):
        voc_metadata_monthly = metadata.iloc[
            metadata.index[metadata["collection_month"] == i].tolist()
        ]
        for index, k in dict(
            collections.Counter(voc_metadata_monthly["lineage"])
        ).items():
            all_voc_df_count.loc[index][i] = k
            all_voc_df_percent.loc[index][i] = round(
                (k / len(metadata.index[metadata["collection_month"] == i].tolist()))
                * 100,
                2,
            )
            all_voc_df_combined.loc[index][
                i
            ] = f"{all_voc_df_count.loc[index][i]} ({all_voc_df_percent.loc[index][i]}%)"
    all_voc_df_count["Total"] = all_voc_df_count.sum(axis=1)
    all_voc_df_count["Percent"] = round(
        (all_voc_df_count["Total"] / all_voc_df_count["Total"].sum()) * 100, 2
    )

    all_voc_df_count.to_excel(writer, f"count")
    all_voc_df_percent.to_excel(writer, f"percent")
    all_voc_df_combined.to_excel(writer, f"combined")
    writer.save()
