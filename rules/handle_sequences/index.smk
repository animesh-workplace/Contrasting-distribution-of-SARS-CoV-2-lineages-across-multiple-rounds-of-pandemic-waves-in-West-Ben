rule handle_sequences:
	message: "Handling sequences for each state"
	input:
		update_log = rules.update.log,
		reference = 'workflow/resources/data/reference.fasta',
	output: directory("{base_path}/General/")
	log: "{base_path}/log/handle_sequences.log"
	threads: 30
	run:
		try:
			metro_states = os.listdir(f"{wildcards.base_path}/Raw_Data")
			for state in metro_states:
				print(f"Fixing metadata of {state}")
				for path, dirs, files in os.walk(f"{wildcards.base_path}/Raw_Data/{state}"):
					if(files):
						for i in files:
							file_url = os.path.join(path, i)
							aligned_url = f"{wildcards.base_path}/General/{state}/alignment/{i.split('.')[0]}_aligned.fasta"
							clade_report_url = f"{wildcards.base_path}/General/{state}/{state}_clade_report.tsv"
							clade_other_url = f"{wildcards.base_path}/General/{state}/other"
							lineage_report_url = f"{wildcards.base_path}/General/{state}/{state}_lineage_report.csv"
							file_type = i.split('_')[0]
							if(file_type == 'sequence'):
								shell(
									f"""
										mkdir -p {wildcards.base_path}/General/{state}/
										cp {file_url} {wildcards.base_path}/General/{state}/
										echo 'Aligning {state} sequences'
										snakemake --snakefile workflow/scripts/aligner/Snakefile --cores {threads} --config path='{wildcards.base_path}/General/{state}/{i.split('.')[0]}'
										echo 'Creating Clade report for {state} sequences'
										nextclade --input-fasta {aligned_url} --output-tsv {clade_report_url} --jobs {threads} \
													--input-root-seq workflow/resources/data/reference.fasta --input-tree workflow/resources/data/tree.json \
													--input-qc-config workflow/resources/data/qc.json --input-pcr-primers workflow/resources/data/primers.csv \
													--input-gene-map workflow/resources/data/genemap.gff --output-dir {clade_other_url}
										echo 'Creating Lineage report for {state} sequences'
										pangolin {aligned_url} --outfile {lineage_report_url} -t {threads}
									"""
								)
		except:
			error_traceback = traceback.format_exc()
			pathlib.Path(str(log)).write_text(error_traceback)
			raise
