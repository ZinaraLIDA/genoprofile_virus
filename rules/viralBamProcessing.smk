rule processViralBam :
	input :
		reportOutput = rules.insertionReport.output,
		fa = rules.virusBreakend.output.fa,
		fai = rules.virusBreakend.output.fai,
		viralBam = rules.virusBreakend.output.viralBam,
		viralBamBai = rules.virusBreakend.output.viralBamBai
	
	params :
		prefix = f"{{prefix}}_chr_{{index}}",
		clippedBam = f"{outputPath}/processViralBam/{{prefix}}_chr_{{index}}_viral_clipped.bam",
		clippedBamBai = f"{outputPath}/processViralBam/{{prefix}}_chr_{{index}}_viral_clipped.bam.bai"
	
	log :
		f"{outputPath}/pipelineLog/processViralBam/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/processViralBam/{{prefix}}_chr_{{index}}.done"
	
	shell :
		"""
		outputDirectory=$(dirname {input.reportOutput})
		#files=$(ls "$outputDirectory"/{params.prefix}*.tsv)
		files=$(find "$outputDirectory" -type f -name {params.prefix}*.tsv)
		if [ -s {input.viralBam} ] && [ -n "$files" ];
			then
				samtools view -h {input.viralBam} | samclip --invert --ref {input.fa} | samtools sort > {params.clippedBam}
				samtools index {params.clippedBam} {params.clippedBamBai}
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""