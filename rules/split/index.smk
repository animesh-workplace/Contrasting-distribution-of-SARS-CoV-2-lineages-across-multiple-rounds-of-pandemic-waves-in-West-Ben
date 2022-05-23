rule split_sequences_metadata:
	message: "Handling metadata and sequences for the wave types"
	input:
		handle_metadata = rules.handle_metadata.output,
		wave_type_url = "workflow/resources/wave_type_date.json"
	output: directory("{base_path}/Output")
	log: "{base_path}/log/output.log"
	threads: 30
	run:
		try:
			with open(input.wave_type_url) as f:
				wave_type = json.loads(f.read())
			metro_states = os.listdir(f"{wildcards.base_path}/Complete/")
			for state in metro_states:
				for path, dirs, files in os.walk(f"{wildcards.base_path}/Complete/{state}"):
					print(f"\nSplitting metadata of {state}")
					metadata_url = f"{wildcards.base_path}/Complete/{state}/complete_metadata_{state}.tsv"
					sequence_url = f"{wildcards.base_path}/Complete/{state}/complete_sequences_{state}.fasta"
					metadata = pandas.read_csv(metadata_url, delimiter = '\t', encoding = 'utf-8', low_memory = False)
					sequence = SeqIO.to_dict(SeqIO.parse(str(sequence_url), 'fasta'))

					shell(
						f"""
							mkdir -p {wildcards.base_path}/Output/{state}/QC_Passed
						"""
					)
					print('\n------------------')		
					print('For QC passed data')
					print('------------------')		
					qc_passed_metadata = metadata.loc[metadata.index[metadata['totalMissing'] < 2990]]
					qc_passed_metadata.reset_index(drop = True, inplace = True)
					qc_passed_metadata = qc_passed_metadata.loc[qc_passed_metadata.index[qc_passed_metadata['location'] != 'Airport']]
					qc_passed_metadata.reset_index(drop = True, inplace = True)
					qc_passed_metadata['date'] = pandas.to_datetime(qc_passed_metadata['date'], format="%Y-%m-%d")
					qc_passed_sequences = [sequence[i] for i in qc_passed_metadata['strain']]
					qc_passed_metadata.to_csv(f"{wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}.tsv", sep = '\t', index = False)
					SeqIO.write(qc_passed_sequences, f"{wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_sequences_{state}.fasta", 'fasta')
					shell(
						f"""
							echo 'Nucleotide Substitution'
							python workflow/scripts/nucleotide_substitution_estimate.py --metadata {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_{state}_nucleotide_substitution_estimate.xlsx
							echo 'Amino acid Substitution'
							python workflow/scripts/amino_acid_substitution_estimate.py --metadata {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_{state}_amino_acid_substitution_estimate.xlsx
							echo 'Lineage Month Report'
							python workflow/scripts/all_lineage_month_report.py --metadata {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}_lineage_report.xlsx
							echo 'Mutation Report'
							python workflow/scripts/mutation_report.py --metadata {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}_mutation_report.xlsx
							echo 'Lineage Substitution Deletion Report'
							python workflow/scripts/lineage_substitution_deletion_report.py --metadata {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/QC_Passed/qc_passed_metadata_{state}_lineage_substitution_deletion_report.tsv
						"""
					)
					for waves in wave_type[state]:
						wave_key = list(waves.keys())[0]
						shell(
							f"""
								mkdir -p {wildcards.base_path}/Output/{state}/{wave_key}
							"""
						)
						print('\n---------------------------')							
						print(f'For {wave_key} passed data')
						print('---------------------------')							
						if(wave_key == 'Wave2'):
							wave_metadata = qc_passed_metadata.loc[(qc_passed_metadata['date'] >= pendulum.parse(waves[wave_key][0]).to_datetime_string())]
						elif(wave_key == 'Pre-Wave1'):
							print(waves[wave_key][1])
							wave_metadata = qc_passed_metadata.loc[(qc_passed_metadata['date'] <= pendulum.parse(waves[wave_key][1]).to_datetime_string())]
						else:
							wave_metadata = qc_passed_metadata.loc[(qc_passed_metadata['date'] >= pendulum.parse(waves[wave_key][0]).to_datetime_string()) & (qc_passed_metadata['date'] <= pendulum.parse(waves[wave_key][1]).to_datetime_string())]
						wave_sequences = [sequence[i] for i in qc_passed_metadata['strain']]
						SeqIO.write(wave_sequences, f"{wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_sequences_{state}.fasta", 'fasta')
						wave_metadata.to_csv(f"{wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_metadata_{state}.tsv", sep = '\t', index = False)
						shell(
							f"""
								echo 'Nucleotide Substitution'
								python workflow/scripts/nucleotide_substitution_estimate.py --metadata {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_{state}_nucleotide_substitution_estimate.xlsx
								echo 'Amino acid Substitution'
								python workflow/scripts/amino_acid_substitution_estimate.py --metadata {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_{state}_amino_acid_substitution_estimate.xlsx
								echo 'Lineage Month Report'
								python workflow/scripts/all_lineage_month_report.py --metadata {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_{state}_lineage_month_report.xlsx
								echo 'Mutation Report'
								python workflow/scripts/mutation_report.py --metadata {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_{state}_mutation_report.xlsx
								echo 'Lineage Substitution Deletion Report'
								python workflow/scripts/lineage_substitution_deletion_report.py --metadata {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_metadata_{state}.tsv --output {wildcards.base_path}/Output/{state}/{wave_key}/{wave_key}_{state}_lineage_substitution_deletion_report.tsv
							"""
						)
		except:
			error_traceback = traceback.format_exc()
			pathlib.Path(str(log)).write_text(error_traceback)
			raise
