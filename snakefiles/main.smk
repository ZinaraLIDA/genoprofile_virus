rulePath=config["RULE_PATH"]
outputPath = config["general_informations"]["OUTPUT"]
workdir : config["general_informations"]["OUTPUT"]

#
reference = config["general_informations"]["FASTA_FILE"]
prefix = config["general_informations"]["PREFIX"]
index = config["general_informations"]["INDEX"]
data = config["general_informations"]["DATA"]

# Import rules
include : rulePath+"/cram2bam.smk"
include : rulePath+"/viralIntegrationDetection.smk"
include : rulePath+"/compressAndIndexVcf.smk"
include : rulePath+"/filterSummary.smk"
include : rulePath+"/filterVcf.smk"
include : rulePath+"/fragmentCovBam.smk"
include : rulePath+"/fragmentCovBed.smk"
include : rulePath+"/getClosestGenes.smk"
include : rulePath+"/insertionReport.smk"
include : rulePath+"/viralBamProcessing.smk"
include : rulePath+"/resultsVisualization.smk"
include : rulePath+"/joinFile.smk"

rule all :
	input :
		expand("{outputPath}/joinedFiles/AllFilesJoined.done", outputPath=outputPath)