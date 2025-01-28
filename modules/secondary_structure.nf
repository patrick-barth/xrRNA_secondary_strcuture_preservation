process predict_secondary_structure {
	tag {mRNA.simpleName}

	input:
	path(mRNA)

	output:
	path("${mRNA.simpleName}.secondary_structure.fold"),	emit: structure,		optional:true
	tuple val(mRNA.simpleName), path("*ss.ps"),				emit: visualization,	optional:true
	path("${task.process}.version.txt"),					emit: version

	"""
	RNAfold --noconv \\
		< ${mRNA} \\
		> ${mRNA.simpleName}.secondary_structure.fold

	echo -e "${task.process}\tRNAfold\t\$(RNAfold --version | cut -d ' ' -f 2)" > ${task.process}.version.txt
	"""
}

process extract_matching_structures {
	tag {sequences.simpleName}
	publishDir "${params.output_dir}", mode: 'copy', pattern: "matching_secondary_structures.fna"

	input:
	path(sequences)

	output:
	path("matching_secondary_structures.fna"),					emit: sequences,	optional:true
	path("${task.process}.version.txt"),						emit: version

	"""
	extract_matching_structures.py --sequence ${sequences} \\
		--output matching_secondary_structures.fna

	echo -e "${task.process}\textract_matching_structures.py\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tpython\t\$(python --version | cut -d ' ' -f 2 )" >> ${task.process}.version.txt
	"""
}

process sort_table {
	tag {query.simpleName}
	publishDir "${params.output_dir}/${query.simpleName}", mode: 'copy', pattern: "${query.simpleName}.sorted.tsv"

	input:
	path(query)

	output:
	path("${query.simpleName}.sorted.tsv"),								emit: sorted_table
	tuple val(query.simpleName), path("${query.simpleName}.sorted.tsv"),	emit: sorted_table_to_merge 
	path("${task.process}.version.txt"),								emit: version

	"""
	sort_table.R \
		--input_file ${query} \
		--output_file ${query.simpleName}.sorted.tsv

	echo -e "${task.process}\sort_table.R\tcustom_script" > ${task.process}.version.txt
    echo -e "${task.process}\tR\t\$(R --version | head -1 | cut -d' ' -f3)" >> ${task.process}.version.txt
	"""
}

process extract_visualizations {
	tag {name}
	publishDir "${params.output_dir}/${name}", mode: 'copy', pattern: "visualization/*ss.ps"

	input:
	tuple val(name), path(query), path(pics)
	val(number_top_hits)

	output:
	path("visualization/*ss.ps")

	"""
	mkdir -p visualization
	tail -n +2 ${query} | head -${number_top_hits} | cut -d\$'\\t' -f1 | tr -d '"' | tr ':' '_' | tr '|' '_' | awk '{print \$0 "_ss.ps"}' > files.txt

	while IFS= read -r file; do
		cp "\$file" visualization/
	done < files.txt
	"""
}