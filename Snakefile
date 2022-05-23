import os
import json
import numpy
import arrow
import base64
import pandas
import pathlib
import fuzzyset
import pendulum
import itertools
import traceback
import subprocess
import collections
from tqdm import tqdm
from Bio import SeqIO
from time import sleep
from Bio.Seq import Seq
from dotenv import load_dotenv
from pandas import ExcelWriter

rule all:
	input:
		uploaded_files = expand(
			"{base_path}/Output",
				base_path=config["base_path"],
		)

include: "rules/update/index.smk"
include: "rules/handle_sequences/index.smk"
include: "rules/handle_metadata/index.smk"
include: "rules/split/index.smk"