import glob

def get_matching_files(wildcards):
	#outputPath = wildcards.outputPath
	prefix = wildcards.prefix
	index = wildcards.index
	files = glob.glob(f"{outputPath}/fragmentCovBam/{prefix}_chr_{index}*.bam")
	return files

rule fragmentCovBed :
	input :
		samtoolsDone = rules.fragmentCovBam.output,
		viralBam = rules.virusBreakend.output.viralBam

	params :
		coordList = f"{outputPath}/filterVcf/{{prefix}}_chr_{{index}}_coord_list.txt",
		samtoolsOutputFile = get_matching_files
	
	log :
		f"{outputPath}/pipelineLog/fragmentCovBed/{{prefix}}_chr_{{index}}.log"

	output :
		f"{outputPath}/fragmentCovBed/{{prefix}}_chr_{{index}}.done"

	shell :
		"""
		outputDirectory=$(dirname {output})
		if [ -s {input.viralBam} ] && [ -s {params.coordList} ];
			then
				for file in {params.samtoolsOutputFile}
					do
						fragment_cov_bed.py \
						--coord_list {params.coordList} \
						--samtools_output_file $file \
						--output_directory $outputDirectory
						touch {output}
					done		
		else
			touch {output}
		fi >> {log} 2>&1
		"""