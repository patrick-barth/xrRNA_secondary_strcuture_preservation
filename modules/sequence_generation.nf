process assemble_sequences {

	input:
	tuple path(xrRNA),path(spacer),path(IRES),path(CDS)
	val(rnd_spacers)
	val(specific_spacers)
	val(percent_specific_spacers)
	val(spacer_min_length)
	val(spacer_max_length)
	val(spacer_fixed_length)
	val(number_sequences)

	output:
	path("generated_sequences.fna"),		emit: sequences
	path("${task.process}.version.txt"),	emit: version

	"""
	generate_sequences.py --xrRNA_file ${xrRNA} \
		--spacer_file ${spacer} \
		--IRES_file ${IRES}	\
		--CDS_file ${CDS} \
		--only_random_spacers ${rnd_spacers} \
		--only_specific_spacers ${specific_spacers} \
		--percent_specific_spacers ${percent_specific_spacers} \
		--spacer_min_len ${spacer_min_len} \
		--spacer_max_len ${spacer_max_len} \
		--spacer_fixed_length ${spacer_fixed_length} \
		--output generated_sequences.fna	\
		--number_sequences ${number_sequences}

	echo -e "${task.process}\tgenerate_sequences.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | cut -d ' ' -f 2 )" >> ${task.process}.version.txt
	"""
}