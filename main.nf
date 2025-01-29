#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//import processes revolving around secondary structure prediction
include{
    predict_secondary_structure
	extract_matching_structures
} from './modules/secondary_structure.nf'

//import processes revolving around sequence manipulation
include{
	assemble_sequences
} from './modules/sequence_generation.nf'

//import general processes
include{
    collect_workflow_metrics
    collect_versions
	collect_files
} from './modules/general_processes.nf'


log.info """\
        ${params.manifest.name} v${params.manifest.version}
        ==========================
        xrRNA file					: ${params.xrRNA_file}
		spacer file					: ${params.spacer_file}
		IRES file					: ${params.ires_file}
		CDS file					: ${params.cds_file}
		----------------------------------------------------------------
		output directory			: ${params.output_dir}
		----------------------------------------------------------------
		Only random spacers			: ${params.only_random_spacers}
		Only specific spacers		: ${params.only_specific_spacers}
		Percent specific spacers	: ${params.percent_specific_spacers}
		Spacer min length			: ${params.spacer_min_length}
		Spacer max length			: ${params.spacer_max_length}
		Spacer fixed length			: ${params.spacer_fixed_length}
		----------------------------------------------------------------
		Number of output sequences	: ${params.number_sequences}
		"""
		.stripIndent()

//input files
input_xrRNA		= Channel.fromPath(params.xrRNA_file)
input_spacer	= Channel.fromPath(params.spacer_file)
input_ires		= Channel.fromPath(params.ires_file)
input_cds		= Channel.fromPath(params.cds_file)

workflow generate_RNA_sequences {
	take:
		input_xrRNA
		input_spacer
		input_ires
		input_cds
	main:
		assemble_sequences(input_xrRNA
			.combine(input_spacer)
			.combine(input_ires)
			.combine(input_cds),
			params.only_random_spacers,
			params.only_specific_spacers,
			params.percent_specific_spacers,
			params.spacer_min_length,
			params.spacer_max_length,
			params.spacer_fixed_length,
			params.number_sequences
		)

		// collect versions
		versions = assemble_sequences.out.version.first()

	emit:
		sequences 	= assemble_sequences.out.sequences
		versions 	= versions
}

workflow test_secondary_structure {
	take:
		rna_sequences
	main:
		predict_secondary_structure(rna_sequences
			.splitFasta(by:10,file:true))
		extract_matching_structures(predict_secondary_structure.out.structure)
		collect_files(extract_matching_structures.out.sequences
			.collect())
		'''sort_table(tabelize_secondary_structures.out.table)
		extract_visualizations(sort_table.out.sorted_table_to_merge
			.join(predict_secondary_structure.out.visualization),
			params.number_top_hits)'''
	
		versions = predict_secondary_structure.out.version.first()
			.concat(extract_matching_structures.out.version.first())

	emit:
		versions = versions
}

workflow {
	generate_RNA_sequences(input_xrRNA,
		input_spacer,
		input_ires,
		input_cds)
	test_secondary_structure(generate_RNA_sequences.out.sequences)


	collect_workflow_metrics()

	// Collect and output tool versions
	versions = generate_RNA_sequences.out.versions
		.concat(test_secondary_structure.out.versions)
	collect_versions(versions.flatten()
		.toList())
}


workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}
