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

metadata = pandas.read_csv(metadata_url, delimiter="\t", encoding="utf-8")
writer = ExcelWriter(output_url)

temp = []
for i in metadata["substitutions"]:
    if isinstance(i, str):
        temp.append(i.split(","))
substitution = list(itertools.chain(*temp))

nucleotide_change = []
for i in substitution:
    nucleotide_change.append(f"{i[0]}->{i[-1]}")

all_changes = dict(collections.Counter(nucleotide_change))

pandas.DataFrame.from_dict(data=[all_changes]).transpose().rename(
    columns={0: "Count"}
).to_excel(writer, "count", header=True)
final_form = pandas.DataFrame(columns=["A", "T", "G", "C"], index=["A", "T", "G", "C"])

# reference_value = {
# 	'A': 8954,
# 	'T': 9594,
# 	'G': 5863,
# 	'C': 5492
# }

total = len(metadata)

for key, value in all_changes.items():
    final_form[key.split("->")[0]][key.split("->")[1]] = value

final_form.to_excel(writer, "final", header=True)
writer.save()
