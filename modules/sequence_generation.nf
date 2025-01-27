process assemble_sequences {

	input:
	tuple path(xrRNA),path(spacer),path(IRES),path(CDS)
	val(rnd_spacer)
	val(specific_spacer)
	val(percent_specific_spacers)
	val(spacer_min_length)
	val(spacer_max_length)
	val(spacer_fixed_length)
	val(number_sequences)

	output:
	path("generated_sequences.fna"),		emit: sequences
	path("${task.process}.version.txt"),	emit: version

	script:
	def spacer_fixed_length = spacer_fixed_length ? '--spacer_fixed_length ' + spacer_fixed_length : ''
	def only_random_spacer	= rnd_spacer ? '--only_random_spacer' : ''
	def only_specific_spacer = specific_spacer ? '--only_specific_spacer' : ''

	"""
	generate_sequences.py --xrRNA_file ${xrRNA} \
		--spacer_file ${spacer} \
		--IRES_file ${IRES}	\
		--CDS_file ${CDS} \
		${only_random_spacer} \
		${only_specific_spacer} \
		--percent_specific_spacers ${percent_specific_spacers} \
		--spacer_min_len ${spacer_min_length} \
		--spacer_max_len ${spacer_max_length} \
		${spacer_fixed_length} \
		--output generated_sequences.fna	\
		--number_sequences ${number_sequences}

	echo -e "${task.process}\tgenerate_sequences.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | cut -d ' ' -f 2 )" >> ${task.process}.version.txt
	"""
}
