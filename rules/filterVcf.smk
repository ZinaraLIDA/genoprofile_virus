rule filterVcf :
	input :
		rules.compressAndIndexVcf.output,
		filterSummaryOutput = rules.filterSummary.output,
		viralBam = rules.virusBreakend.output.viralBam,
		vcf = rules.virusBreakend.output.vcf

	params :
		filteredSummary = f"{outputPath}/filterSummary/{{prefix}}_chr_{{index}}_markdup_all_virus.tsv"
	
	log :
		f"{outputPath}/pipelineLog/filterVcf/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/filterVcf/{{prefix}}_chr_{{index}}.done"

	shell :
		"""
		outputDirectory=$(dirname {output})
		
		if [ -s {input.viralBam} ] && [ -f {params.filteredSummary} ];
			then
				filter_vcf.py \
				--vcf {input.vcf}.gz \
				--filtered_summary {params.filteredSummary} \
				--output_directory $outputDirectory
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""