rule cram2bam :
	input :
		cram = f"{data}/{{prefix}}_chr_{{index}}_markdup.cram",
		crai = f"{data}/{{prefix}}_chr_{{index}}_markdup.cram.crai"
	
	params :
		reference = config["general_informations"]["FASTA_FILE"]
	
	log :
		f"{outputPath}/pipelineLog/bam2cram/{{prefix}}_chr_{{index}}.log"
	
	output :
		bam = f"{outputPath}/bam/{{prefix}}_chr_{{index}}_markdup.bam",
		bai = f"{outputPath}/bam/{{prefix}}_chr_{{index}}_markdup.bam.bai"
	
	shell :
		"""
		samtools view -b -T {params.reference} -o {output.bam} {input.cram} >> {log} 2>&1
		samtools index -b {output.bam} >> {log} 2>&1
		"""