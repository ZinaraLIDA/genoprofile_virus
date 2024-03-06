rule insertionReport :
	input :
		rules.getClosestGenes.output,
		viralBam = rules.virusBreakend.output.viralBam
	
	params :
		UpdatedCoordList = f"{outputPath}/fragmentCovBed/{{prefix}}_chr_{{index}}_updated_coord_list.txt",
		closestTsv = f"{outputPath}/getClosestGenes/{{prefix}}_chr_{{index}}.gene_distance.tsv"

	log :
		f"{outputPath}/pipelineLog/insertionReport/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/insertionReport/{{prefix}}_chr_{{index}}.done"

	shell :
		"""
		outputDirectory=$(dirname {output})
		if [ -s {input.viralBam} ] && [ -f {params.UpdatedCoordList} ];
			then
				insertion_report.py \
				--updated_coord_list {params.UpdatedCoordList} \
				--closest_tsv {params.closestTsv}  \
				--output_directory $outputDirectory
				touch {output}
		else
			touch {output}
		fi >> {log} 2>&1
		"""