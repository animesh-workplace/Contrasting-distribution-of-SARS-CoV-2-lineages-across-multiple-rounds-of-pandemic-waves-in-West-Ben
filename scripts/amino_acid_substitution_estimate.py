import pandas
import argparse
import itertools
import collections
from tqdm import tqdm
from pandas import ExcelWriter

parser = argparse.ArgumentParser()
parser.add_argument("--metadata", help="Enter Metadata")
parser.add_argument("--output", help="Enter output file")
args = parser.parse_args()
metadata_url = args.metadata
output_url = args.output

genes = [
    "E",
    "M",
    "N",
    "ORF1a",
    "ORF1b",
    "ORF3a",
    "ORF6",
    "ORF7a",
    "ORF7b",
    "ORF8",
    "ORF9b",
    "S",
]
metadata = pandas.read_csv(metadata_url, delimiter="\t", encoding="utf-8")
writer = ExcelWriter(output_url)

temp = []
for i in metadata["aaSubstitutions"]:
    if isinstance(i, str):
        temp.append(i.split(","))
aaSubstitution = list(itertools.chain(*temp))

final_data = {}
for i in genes:
    aaSubstitution_pos = []
    for j in aaSubstitution:
        if j.startswith(f"{i}"):
            aaSubstitution_pos.append(j.split(":")[1][1:-1])
    substitution_pos_count = dict(collections.Counter(aaSubstitution_pos))
    pandas.DataFrame.from_dict(data=[substitution_pos_count]).transpose().rename(
        columns={0: "Count"}
    ).to_excel(writer, i, header=True)
    final_data[i] = len(substitution_pos_count.keys())

pandas.DataFrame.from_dict(data=[final_data]).transpose().rename(
    columns={0: "Count"}
).to_excel(writer, "total", header=True)
writer.save()
