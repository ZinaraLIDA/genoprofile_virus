rule compressAndIndexVcf :
	input :
		vcf = rules.virusBreakend.output.vcf,
		viralBam = rules.virusBreakend.output.viralBam
	
	log :
		f"{outputPath}/pipelineLog/compressAndIndexVcf/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/compressAndIndexVcf/{{prefix}}_chr_{{index}}.done"

	shell :
		"""
		if [ -s {input.viralBam} ];
			then
				bgzip -fc {input.vcf} > {input.vcf}.gz
				tabix -p vcf {input.vcf}.gz
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""