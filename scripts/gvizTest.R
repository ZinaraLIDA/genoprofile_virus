library(GenomicRanges)
library("Gviz")
library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)
library(VariantAnnotation)
library(lattice)
library("optparse")
library("seqinr")

write("libary loaded", stdout())

#Fix grid plot for when many insersion occur

 
#Input arguments
option_list = list(
	make_option(c("--fullbam"), type="character", default=NULL, 
				help="bam of aligned reads against viral ref", metavar="character"),
	make_option(c("--clipbam"), type="character", default=NULL, 
				help="bam of soft clipped reads aligned against viral ref", metavar="character"),
	make_option(c("--fastaref"), type="character", default=NULL, 
				help="fasta containing the adjusted viral sequence made by virus breakend", metavar="character"),
	make_option(c("--vcf"), type="character", default=NULL, 
				help="virus breakend result vcf file", metavar="character"),
	make_option(c("--filteredres"), type="character", default=NULL, 
				help="filtered virus table from virusBreakend", metavar="character"),
	make_option(c("--sample"), type="character", default=NULL, 
				help="sample name", metavar="character"),
	make_option(c("--o"), type="character", default=NULL, 
				help="output directory", metavar="character"),
	make_option(c("--bands"), type="character", default=NULL, 
				help="Location of each cytoBands of Homo Sapiens genome build hg38", metavar="character")

#__________________________________________________________
	# make_option(c("--gffDir"), type="character", default=NULL, 
	#			 help="directory for the viral gff annotation files", metavar="character")
#__________________________________________________________
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sample = opt$sample

outdir <- opt$o

#General
# Disable force conventional naming for chromosome (e.g. "chr3")
options(ucscChromosomeNames=FALSE)


#Infput files
#BAM
# Reads aligned to viral reference
full_bam <- opt$fullbam
# Soft clipped read (visualize read that support integration)
clipped_bam <- opt$clipbam

#__________________________________________________________
#GFF
#gffDir <- opt$gffDir
#__________________________________________________________


#Fasta
#Fasta for the viral reference
fastaRef <- opt$fastaref
#VCF
#Vcf file with insersion coordinates in both host and viral reference
vcfFile <- opt$vcf
vcf_df <- readVcf(vcfFile)

#Chromosome
# Name of the viral reference in the vcf, fastas and GFF
#import filter result
res_table <- read.table(file = opt$filteredres, sep = '\t', header = TRUE)

#get chromosome names (viral sequence here)
chr_list <- list(res_table$rname)

# for each chromosome (viral sequence here) plot insertion

#filter virus with insertion detected
res_table <- filter(res_table, integrations > 0)

#
cytoBands = read.table(opt$bands, header = TRUE, sep = "\t")

for (row in 1:nrow(res_table)) {
	chr <- res_table[row, 'rname']
	gffname <- gsub("adjusted_", "", res_table$rname)
	vir_name <- res_table[row, 'name_assigned']
	vir_name <- gsub(" ", "_", vir_name)
	cat("check")
	
	#Input Bam
	#All
	peakReads <- AlignmentsTrack(full_bam, isPaired=TRUE, name='coverage')
	#Soft clip (to modify)
	sec_reads <- AlignmentsTrack(clipped_bam, isPaired=TRUE, name='split read')

	# Reference track
	#Parse GFF
	# Parse viral reference too extract start end and id of the viral genes

	#Import gff
	#gffFile <- paste(gffDir, '/', gffname, '.gff' )
	#___________________________________________________________________________
	# df <- readGFF(gffFile)
	# df1 <- as_tibble(df)
	# df2 <- df1 %>% filter(type == 'gene')
	# start_ref<-pull(df2, start)
	# end_ref<-pull(df2,end)
	# id_ref<-pull(df2,Name)
	# strand_ref<-pull(df2,strand)
	# aTrack <- AnnotationTrack(start=start_ref, end=end_ref, strand = strand_ref, chromosome=chr, id=id_ref, name='Annotation', stacking = "full")
	#___________________________________________________________________________
	
	#Parse VCF
	#Parse VCF to extract all the breakpoints in the viral genome
	vcf_df1 <- as_tibble(rowRanges(vcf_df))
	if (length(vcf_df1) == 10) {
		print(sprintf("vcf file of %s is not empty", sample))
		#print(length(vcf_df1))
		vcf_df2 <- vcf_df1 %>% filter(seqnames == chr)
		bp_start <- pull(vcf_df2, start)
		bp_width <- pull(vcf_df2, width)
		#Combine the 2 alignment tracks to display the viral insertion location on both (vertical red lines)
		ht <- HighlightTrack(trackList = list(peakReads, sec_reads), start = bp_start, width = bp_width, chromosome = chr, inBackground=FALSE)
		
		#Genome axis track to display coordinates
		genomeAxis <- GenomeAxisTrack(name="MyAxis")

		#Ref sequence track
		fasta_obj <- read.fasta(file = fastaRef, as.string = T)
		ref_start <- 0
		ref_end <- getLength(fasta_obj[chr])


		#Plot without annotation
		write("TYPE:", stdout())
		write(typeof(sample), stdout())
		out_f_fname <- paste(outdir, '/', sample, '_', vir_name, '.png', sep="")
		print("ok6")
		write(out_f_fname, stdout())
		png(file=out_f_fname,
		width=600, height=350)
		#__________________________________________________________
		plotTracks(list(genomeAxis, ht), from=ref_start, to=ref_end, chromosome=chr, featureAnnotation = "id",sizes=c(1,4,4))
		#__________________________________________________________
		#Change when gff is added (and remove line above)
		#__________________________________________________________
		#plotTracks(list(genomeAxis, aTrack, ht), from=ref_start, to=ref_end, chromosome=chr, type='coverage', featureAnnotation = "id",sizes=c(1,1,4,4))
		#__________________________________________________________
		dev.off()


		#Human chromosome view
		#filter vcf_df
		#Parse vcf to get insertion in the human genome
		chrom_list = seq(1,23,by=1)
		append(chrom_list, c('x', 'y'))
		vcf_df3 <- vcf_df1 %>% filter(seqnames %in% chrom_list) %>% mutate(seqnames = paste("chr",seqnames, sep=""))
		write(typeof(sample), stdout())
		out_f_fname <- paste(outdir, '/', sample, '_', vir_name, '_human', '.png', sep="")
		write(out_f_fname, stdout())
		png(file=out_f_fname,
		width=600, height=350)
		#grid plot for all insersion in human chromosome
		ncols <- 1
		nrows <- nrow(vcf_df3)
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrows, ncols)))
		for(i in 1:nrows) {
		pushViewport(viewport(layout.pos.col = 1,
								layout.pos.row = i))
		print(vcf_df3)
		print(vcf_df3[, 1])
		itrack <- IdeogramTrack(
			genome = "GRCh38.92", chromosome = vcf_df3[i,1],
			bands = cytoBands
		)
		plotTracks(itrack, from=vcf_df3[i,2], to=vcf_df3[i,3], add = TRUE )
		popViewport(1)
		}
		dev.off()
	} else {
		print(sprintf("vcf file of %s is empty", sample))
	}
}






	# TODO: reimplement when viral annotation is added
	#Plot
	#Display reads or coverage?
	#with ref to ad SNPs and consensus in the coverage track
	#plotTracks(list(genomeAxis, aTrack, ht, strack), from=ref_start, to=ref_end, chromosome=chr, type='coverage', featureAnnotation = "id",sizes=c(1,1,4,4,1))
	#without ref
	#plotTracks(list(genomeAxis, aTrack, ht), from=ref_start, to=ref_end, chromosome=chr, type='coverage', featureAnnotation = "id",sizes=c(1,1,4,4))
	#Read view
	#plotTracks(list(genomeAxis, aTrack, ht), from=ref_start, to=ref_end, chromosome=chr, featureAnnotation = "id",sizes=c(1,1,4,4))