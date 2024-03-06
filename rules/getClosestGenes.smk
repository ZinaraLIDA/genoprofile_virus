rule getClosestGenes :
	input :
		viralBam = rules.virusBreakend.output.viralBam,
		OutputFragmentCovBed = rules.fragmentCovBed.output
	
	params :
		UpdatedCoordList = f"{outputPath}/fragmentCovBed/{{prefix}}_chr_{{index}}_updated_coord_list.txt",
		bed = config["getClosestGenes"]["BED"]
	
	log :
		f"{outputPath}/pipelineLog/getClosestGenes/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/getClosestGenes/{{prefix}}_chr_{{index}}.done"

	shell :
		"""
		outputDirectory=$(dirname {output})
		if [ -s {input.viralBam} ] && [ -f {params.UpdatedCoordList} ];
			then
				get_closest_genes.py \
				--updated_coord_list {params.UpdatedCoordList} \
				--bed {params.bed}  \
				--output_directory $outputDirectory
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""