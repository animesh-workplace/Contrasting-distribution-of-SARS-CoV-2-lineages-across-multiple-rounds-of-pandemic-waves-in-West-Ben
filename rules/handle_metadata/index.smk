rule handle_metadata:
	message: "Handling metadata for each state"
	input:
		handle_sequences = rules.handle_sequences.output,
		reference = 'workflow/resources/data/reference.fasta',
		indian_region = 'workflow/resources/indian_region.tsv'
	output: directory("{base_path}/Complete")
	log: "{base_path}/log/handle_metadata.log"
	threads: 30
	run:
		try:
			metro_states = os.listdir(f"{wildcards.base_path}/Raw_Data")
			for state in metro_states:
				print(f"Combining reports and metadata of {state}")
				for path, dirs, files in os.walk(f"{wildcards.base_path}/Raw_Data/{state}"):
					if(files):
						for i in files:
							nextstrain_labels = ['strain', 'virus', 'gisaid_epi_isl', 'genbank_accession', 'date', 'region', 'country', 'division', 'location', 'region_exposure', 'country_exposure', 'division_exposure', 'segment', 'length', 'host', 'age', 'sex', 'originating_lab', 'submitting_lab', 'authors', 'url', 'title', 'paper_url', 'date_submitted', 'purpose_of_sequencing']
							file_url = os.path.join(path, i)
							aligned_url = f"{wildcards.base_path}/General/{state}/alignment/{i.split('.')[0].replace('metadata', 'sequence')}_aligned.fasta"
							clade_report_url = f"{wildcards.base_path}/General/{state}/{state}_clade_report.tsv"
							lineage_report_url = f"{wildcards.base_path}/General/{state}/{state}_lineage_report.csv"
							output_base_url = f"{wildcards.base_path}/Complete/{state}/"
							file_type = i.split('_')[0]
							if(file_type == 'metadata'):
								metadata = pandas.read_csv(file_url, delimiter = '\t', encoding = 'utf-8', low_memory = False)
								sequence = SeqIO.to_dict(SeqIO.parse(str(aligned_url), 'fasta'))

								nextclade_metadata = pandas.read_csv(clade_report_url, delimiter = '\t', encoding = 'utf-8', low_memory = False)
								nextclade_metadata.rename(columns = {'seqName': 'strain'}, inplace = True)

								pangolin_metadata = pandas.read_csv(lineage_report_url, delimiter = ',', encoding = 'utf-8', low_memory = False)
								pangolin_metadata.rename(columns = {'taxon': 'strain'}, inplace = True)

								nextclade_pangolin = pandas.merge(
									nextclade_metadata[['strain', 'clade', 'totalInsertions', 'totalMissing', 'totalNonACGTNs', 'nonACGTNs', 'substitutions', 'deletions', 'aaSubstitutions', 'aaDeletions']],
									pangolin_metadata[['strain', 'lineage', 'scorpio_call', 'note', 'pangoLEARN_version']],
									on = 'strain', how = 'inner'
								)

								region_type = pandas.read_csv(input.indian_region, delimiter = '\t', encoding = 'utf-8').set_index('State').T.to_dict()
								
								# For Nextstrain Analysis
								nextstrain_metadata = pandas.DataFrame(columns = nextstrain_labels)
								nextstrain_metadata = nextstrain_metadata.assign(
									strain = metadata['Virus name'],
									lab_id = metadata['Sample ID given by the submitting lab'],
									last_vaccinated = metadata['Last vaccinated'],
									virus = metadata['Type'],
									gisaid_epi_isl = [f'EPI_ISL_{i}' for i in metadata.index],
									genbank_accession = ['?' for i in metadata.index],
									date = metadata['Collection date'],
									region = metadata['Location'],
									country = metadata['Country'],
									division = metadata['State'],
									location = metadata['District'],
									region_type = [region_type[i]['Region'] for i in metadata['State']],
									region_exposure = metadata['Location'],
									country_exposure = metadata['Country'],
									division_exposure = metadata['State'],
									segment = ['genome' for i in metadata.index],
									length = [len(sequence[i]) for i in metadata['Virus name']],
									host = metadata['Host'],
									age = metadata['Patient age'],
									sex = metadata['Gender'],
									originating_lab = metadata['Originating lab'],
									submitting_lab = metadata['Submitting lab'],
									authors = metadata['Authors'],
									url = ['?' for i in metadata.index],
									title = ['?' for i in metadata.index],
									paper_url = ['?' for i in metadata.index],
									date_submitted = metadata['Collection date'],
									purpose_of_sequencing = ['?' for i in metadata.index]
								)								

								shell(
									f"""
										mkdir -p {output_base_url}
										cp {aligned_url} {wildcards.base_path}/Complete/{state}/complete_sequences_{state}.fasta
									"""
								)

								nextstrain_metadata = nextstrain_metadata.merge(nextclade_pangolin, on = 'strain', how = 'inner')
								nextstrain_metadata.to_csv(f"{output_base_url}/complete_metadata_{state}.tsv", sep = '\t', index = False)
		except:
			error_traceback = traceback.format_exc()
			pathlib.Path(str(log)).write_text(error_traceback)
			raise
