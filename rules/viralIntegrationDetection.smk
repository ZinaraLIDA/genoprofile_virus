rule virusBreakend :
	input :
		f"{outputPath}/bam/{{prefix}}_chr_{{index}}_markdup.bam"
	
	params :
		kraken2Db = config["virusBreakend"]["KRAKEN2DB_PATH"],
		options = config["virusBreakend"]["OPTIONS"],
		reference = config["general_informations"]["FASTA_FILE"],
		jarPath = config["virusBreakend"]["JAR_PATH"]
	
	log :
		f"{outputPath}/pipelineLog/virusbreakend/{{prefix}}_chr_{{index}}.log"
	
	output :
		viralBam = touch(f"{outputPath}/virusbreakend/{{prefix}}_chr_{{index}}_markdup.vcf.virusbreakend.working/adjusted/{{prefix}}_chr_{{index}}_markdup.bam.viral.bam"),
		viralBamBai = touch(f"{outputPath}/virusbreakend/{{prefix}}_chr_{{index}}_markdup.vcf.virusbreakend.working/adjusted/{{prefix}}_chr_{{index}}_markdup.bam.viral.bam.bai"),
		fa = touch(f"{outputPath}/virusbreakend/{{prefix}}_chr_{{index}}_markdup.vcf.virusbreakend.working/adjusted/{{prefix}}_chr_{{index}}_markdup.vcf.viral.fa"),
		fai = touch(f"{outputPath}/virusbreakend/{{prefix}}_chr_{{index}}_markdup.vcf.virusbreakend.working/adjusted/{{prefix}}_chr_{{index}}_markdup.vcf.viral.fa.fai"),
		vcf = f"{outputPath}/virusbreakend/{{prefix}}_chr_{{index}}_markdup.vcf",
		tsv = f"{outputPath}/virusbreakend/{{prefix}}_chr_{{index}}_markdup.vcf.summary.tsv"
	
	shell :
		"""
		BASEDIR=$(dirname {output.vcf})
		virusbreakend -r {params.reference} \
		-o {output.vcf} \
		-w "$BASEDIR" \
		--kraken2db {params.kraken2Db} \
		--jar {params.jarPath} \
		{params.options} \
		{input} >> {log} 2>&1
		"""