rule filterSummary :
	input :
		summary = rules.virusBreakend.output.tsv,
		viralBam = rules.virusBreakend.output.viralBam
	
	params :
		virusList = config["filterSummary"]["VIRUS_LIST"]
	
	log :
		f"{outputPath}/pipelineLog/filterSummary/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/filterSummary/{{prefix}}_chr_{{index}}.done"

	shell :
		"""
		outputDirectory=$(dirname {output})
		if [ -s {input.viralBam} ];
			then
				filter_summary.py \
				--summary {input.summary} \
				--virus_list {params.virusList} \
				--output_directory $outputDirectory
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""