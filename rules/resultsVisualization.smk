rule visualizeResults :
	input :
		rules.processViralBam.output,
		fa = rules.virusBreakend.output.fa,
		fai = rules.virusBreakend.output.fai,
		viralBam = rules.virusBreakend.output.viralBam,
		viralBamBai = rules.virusBreakend.output.viralBamBai,
		vcf = rules.virusBreakend.output.vcf
	
	params :
		bands = config["visualizeResults"]["BANDS"],
		tableAllVirus = f"{outputPath}/filterSummary/{{prefix}}_chr_{{index}}_markdup_all_virus.tsv",
		clippedBam = rules.processViralBam.params.clippedBam,
		clippedBamBai = rules.processViralBam.params.clippedBamBai
	
	log :
		f"{outputPath}/pipelineLog/visualizeResults/{{prefix}}_chr_{{index}}.log"
	
	output :
		f"{outputPath}/visualizeResults/{{prefix}}_chr_{{index}}.done"
	
	shell :
		"""
		outputDirectory=$(dirname {output})
		filename=$(basename {input.vcf})
		sample=$(echo "$filename" | sed 's/\.vcf$//')
		if [ -f {params.clippedBam} ] && [ -f {params.clippedBamBai} ]
			then
				virusbreakendr.R \
				--fullbam {input.viralBam} \
				--clipbam {params.clippedBam} \
				--fastaref {input.fa} \
				--vcf {input.vcf} \
				--filteredres {params.tableAllVirus} \
				--sample $sample \
				--o $outputDirectory \
				--bands {params.bands}
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""