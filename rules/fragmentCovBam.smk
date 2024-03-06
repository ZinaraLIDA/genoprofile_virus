rule fragmentCovBam :
	input :
		rules.filterVcf.output,
		bam = f"{outputPath}/bam/{{prefix}}_chr_{{index}}_markdup.bam",
		viralBam = rules.virusBreakend.output.viralBam

	params :
		coordList = f"{outputPath}/filterVcf/{{prefix}}_chr_{{index}}_coord_list.txt"
	
	log :
		f"{outputPath}/pipelineLog/fragmentCovBam/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/fragmentCovBam/{{prefix}}_chr_{{index}}.done"

	shell :
		"""
		outputDirectory=$(dirname {output})
		if [ -s {input.viralBam} ] && [ -s {params.coordList} ];
			then
				fragment_cov_bam.py \
				--bam {input.bam} \
				--coord_list {params.coordList} \
				--output_directory $outputDirectory
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""