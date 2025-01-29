/*
 * Collect all files from channel and emit a single file
 */
process collect_files{
    publishDir "${params.output_dir}", mode: "move", pattern: "collected_valid_secondary_structures.fna"

    input:
    path(sequences)

    output:
    path("collected_valid_secondary_structures.fna"), emit: collected_files

    """
    cat ${sequences} > collected_valid_secondary_structures.fna
    """
}

/*
 * Collects a variety of different metrics such as the command line used to execute the workflow and the parameters used
 */
process collect_workflow_metrics{
	publishDir "${params.output_dir}/metadata", mode: "move"

	output:
	path("workflow_metrics.txt"), emit: output

	"""
	cat <<EOF > workflow_metrics.txt
    Author: ${params.manifest.author}
    Pipeline version: ${params.manifest.version}
    Nextflow version: ${nextflow.version}
    Working directory: ${workflow.workDir}
    Project directory: ${workflow.projectDir}
    User name: ${workflow.userName}
    Execution directory: ${workflow.launchDir}
    Executed script: ${workflow.scriptFile}
    Configuration file: ${workflow.configFiles}
    Configuration profile: ${workflow.profile}
    Command line: ${workflow.commandLine}
    Parameters: ${params}
    Unique session ID: ${workflow.sessionId}

    Container engine: ${workflow.containerEngine}    
    Containers used: ${workflow.container}

    Git repository: ${workflow.repository}
    Repository revision: ${workflow.revision}
    EOF
	"""
}

/*
 * Collects all version outputs and merges them into a single file
 */
process collect_versions {
    publishDir "${params.output_dir}/metadata", mode: 'copy', pattern: "tool_versions.txt"

    input:
    path(query)

    output:
    path('tool_versions.txt'), emit: txt_to_output_dir

    script:
    """
    echo -e "Process\ttool\tversion" > tool_versions.txt
    for i in ${query}
	do
		cat \$i >> tool_versions.txt
	done
    """
}