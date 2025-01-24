#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//import processes revolving around secondary structure prediction
include{
    predict_secondary_structure
	tabelize_secondary_structures
	sort_table
	extract_visualizations
} from './modules/secondary_structure.nf'

//import processes revolving around sequence manipulation
include{
	split_target_mRNA
} from './modules/sequence_generation.nf'

//import general processes
include{
    collect_workflow_metrics
    collect_versions
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
//TODO: Adapt everything below
input_mRNA		= Channel.fromPath(params.target_mRNA)
input_stem		= Channel.fromPath(params.stem_file)
input_linker	= Channel.fromPath(params.linker_file)

workflow generate_circasRNAs {
	take:
		input_mRNA
		input_stem
		input_linker
	main:
		predict_secondary_structure(input_mRNA)
		split_target_mRNA(predict_secondary_structure.out.structure
			.combine(input_stem)
			.combine(input_linker),
			params.loop_length)

		// collect versions
		versions = predict_secondary_structure.out.version.first()
			.concat(split_target_mRNA.out.version.first())

	emit:
		split_mRNA 	= split_target_mRNA.out.split_seq
		versions = versions
}

workflow test_secondary_structure {
	take:
		split_RNAs
	main:
		predict_secondary_structure(split_RNAs)
		tabelize_secondary_structures(predict_secondary_structure.out.structure)
		sort_table(tabelize_secondary_structures.out.table)
		extract_visualizations(sort_table.out.sorted_table_to_merge
			.join(predict_secondary_structure.out.visualization),
			params.number_top_hits)
	
		versions = predict_secondary_structure.out.version.first()
			.concat(tabelize_secondary_structures.out.version.first())
			.concat(sort_table.out.version.first())
	emit:
		versions = versions
}

workflow {
	generate_circasRNAs(input_mRNA,
		input_stem,
		input_linker)
	test_secondary_structure(generate_circasRNAs.out.split_mRNA)


	collect_workflow_metrics()

	// Collect and output tool versions
	versions = generate_circasRNAs.out.versions
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