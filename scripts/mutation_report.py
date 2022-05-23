import pandas
import argparse
import itertools
import collections
from pandas import ExcelWriter

parser = argparse.ArgumentParser()
parser.add_argument("--metadata", help="Enter Nexstrain metadata")
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

temp = []
for i in metadata["aaDeletions"]:
    if isinstance(i, str):
        temp.append(i.split(","))
aaDeletions = list(itertools.chain(*temp))
top_aaDeletions = dict(collections.Counter(aaDeletions))
df = (
    pandas.DataFrame.from_dict(data=[top_aaDeletions])
    .transpose()
    .rename(columns={0: "Count"})
)
df["percent"] = round((df["Count"] / df["Count"].sum()) * 100, 2)
df.to_excel(writer, "deletions", header=True)
for_paper = df.loc[df.index[df["percent"] > 5]]

temp = []
for i in metadata["substitutions"]:
    if isinstance(i, str):
        temp.append(i.split(","))
aaDeletions = list(itertools.chain(*temp))
top_aaDeletions = dict(collections.Counter(aaDeletions))
df = (
    pandas.DataFrame.from_dict(data=[top_aaDeletions])
    .transpose()
    .rename(columns={0: "Count"})
)
df["percent"] = round((df["Count"] / df["Count"].sum()) * 100, 2)
df.to_excel(writer, "nucSubstitution", header=True)
for_paper = pandas.concat([for_paper, df.loc[df.index[df["percent"] > 5]]])

temp = []
for i in metadata["deletions"]:
    if isinstance(i, str):
        temp.append(i.split(","))
aaDeletions = list(itertools.chain(*temp))
top_aaDeletions = dict(collections.Counter(aaDeletions))
df = (
    pandas.DataFrame.from_dict(data=[top_aaDeletions])
    .transpose()
    .rename(columns={0: "Count"})
)
df["percent"] = round((df["Count"] / df["Count"].sum()) * 100, 2)
df.to_excel(writer, "nucDeletions", header=True)
for_paper = pandas.concat([for_paper, df.loc[df.index[df["percent"] > 5]]])

for i in genes:
    gene_aaSubstitution = []
    for j in aaSubstitution:
        if j.startswith(f"{i}"):
            gene_aaSubstitution.append(j)
    gene_subsitution_counter = dict(collections.Counter(gene_aaSubstitution))
    df = (
        pandas.DataFrame.from_dict(data=[gene_subsitution_counter])
        .transpose()
        .rename(columns={0: "Count"})
    )
    df["percent"] = round((df["Count"] / df["Count"].sum()) * 100, 2)
    df.to_excel(writer, i, header=True)
    for_paper = pandas.concat([for_paper, df.loc[df.index[df["percent"] > 5]]])

all_substitution = dict(collections.Counter(aaSubstitution))
df = (
    pandas.DataFrame.from_dict(data=[all_substitution])
    .transpose()
    .rename(columns={0: "Count"})
)
df["percent"] = round((df["Count"] / df["Count"].sum()) * 100, 2)
df.to_excel(writer, "All", header=True)
for_paper.to_excel(writer, "Paper", header=True)

writer.save()
