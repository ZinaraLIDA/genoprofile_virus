#!/usr/bin/python3

import glob
import argparse
import re
import os
import subprocess
import sys

import pandas as pd
from pysam import VariantFile

#requirments: samtools and bedtools

#TODO add cpus to samtools

parser = argparse.ArgumentParser()

parser.add_argument("-vcf", type=str,
					help="Path to virus breaken output VCF")
parser.add_argument("-bam", type=str,
					help="Path to the star/bwa BAM file")
parser.add_argument("-bed", type=str,
					help="Path to the bed annotation file from the human genome reference")
parser.add_argument("-summary", type=str,
					help="Path to the virus breakend summary tsv output")
parser.add_argument("-virus_list", type=str,
					help="Path to the list for the viruses of intrest")
parser.add_argument("-o", type=str,
					help="Output directory")

args = parser.parse_args()

def filter_summary(summary, virus_list, output_directory):
# Filter tsv for virus of interest and extract relevant fields
	summary_df = pd.read_csv(summary, sep = '\t', header = 0)
	#sample = summary.split("/")[-1].split('.')[0]
	sample = summary.split("/")[-1].replace('.vcf.summary.tsv', '')
	print(sample)
	filtered_df = summary_df[(summary_df["name_species"].isin(virus_list)) | (summary_df["name_genus"].isin(virus_list))][["name_genus", "name_species", "name_assigned", "rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "integrations"]]
	if not filtered_df.empty :
		filtered_df.to_csv(f"{output_directory}/{sample}_all_virus.tsv", sep = '\t', index=False)
	else :
		filtered_df = pd.DataFrame()
		with open(f"{output_directory}/{sample}_all_virus.tsv", "w") as filin :
			pass
	return filtered_df
	
# To be called for each line of filtered tsv for wich intergration was detected
def filter_vcf(vcf, rname):
# Extract BVF and host genome coordinates for each integration
	coord_list = []
	in_vcf = VariantFile(vcf)
	# Get virus name
	for rec in in_vcf.fetch(rname):
		bvf = rec.info["BVF"]
		alt = rec.alts[0]
		pos_re = r"([0-9]+|[X,Y]):([0-9]+)"
		m = re.search(pos_re, alt)
		chrom = m.group(1)
		pos = m.group(2)
		coord_list.append([chrom, pos, bvf])
	#TOD: define best format
	return coord_list

def fragment_cov_bam(bam, coord_list, window_size = 4000):
	sample_name = bam.split("/")[-1].strip(".bam")
	#extract for each insersion
	#use dictionary next time -_-
	updated_coord_list = []
	for coord in coord_list:
		chrom = coord[0]
		pos = int(coord[1])
		bvf = coord[2]
		start = int(pos - (window_size/2))
		end = int(pos + (window_size/2))
		updated_coord = [chrom, pos, bvf]
		#Filter fully matched pairs?
		#compute fragment coverage
		#would output only the coordlist with added frgment coverage
		cmd = f'samtools view -hb {bam} "{chrom}:{start}-{end}" | bedtools genomecov -ibam stdin -pc -bga'
		result = subprocess.check_output(cmd, shell=True)
		cov_df = pd.DataFrame([row.split('\t') for row in result.decode().split('\n')], 
			columns=['chr', 'start', 'end', 'cov'])
		#Last row is only none
		cov_df.dropna(axis="index", inplace = True )
		cov_df = cov_df[:-1].astype({'start':int, 'end':int, 'cov':int})
		df = cov_df[(pos >= cov_df['start']) & (pos < cov_df['end']) & (cov_df['chr'] == chrom)]
		cov = int(df['cov'])
		updated_coord.append(cov)
		updated_coord.append(float(bvf/(bvf+cov)))
		updated_coord_list.append(updated_coord)
	for coord in updated_coord_list:
		coord = [int(x) for x in coord[:-1] if x not in ['X', 'Y']]
	updated_coord_list.sort()
	return updated_coord_list, sample_name

def closest_genes(coord_list, bed, output_directory, sample):
	#uses bedtools closest to find genes around the insersion
	#Create bed file for the insertion
	with open(f'{output_directory}/{sample}.insertion.bed', 'w') as out_bed:
		for elem in coord_list:
			bed_line = elem[:2]
			bed_line.append(int(bed_line[-1])+1)
			bed_line = [str(x) for x in bed_line]
			out_bed.write('\t'.join(bed_line) + '\n')
			print(out_bed)
	#intersect with bed tools
	cmd = f'bedtools closest -a {output_directory}/{sample}.insertion.bed -b {bed} -k 5 -D ref > {output_directory}/{sample}.gene_distance.tsv'
	os.system(cmd)
	print('bedtools done')

def insertion_report(coord_list, sample, virus, output_directory):
	columns = ['ins_chr', 'ins_start', 'ins_start1', 'an_chr', 'an_start', 'an_end' , 'an_id', 'dist']
	res_columns = ['Chr', 'Pos', 'supporting_frag', 'host_frag', 'vaf', 'interupeted_gene', 'nearby_genes']
	closest_tsv = f"{output_directory}/{sample}.gene_distance.tsv"
	distance_df = pd.read_csv(closest_tsv, sep = '\t', header = None, names = columns)
	res_list =[]
	for coord in coord_list:
		sub_df = distance_df[(distance_df["ins_chr"].astype('str') == coord[0]) &
		 (distance_df["ins_start"] == coord[1])]
		gene_list = sub_df['an_id'].to_list()
		distance_list = sub_df['dist'].to_list()
		genes_str =''
		for i in range(len(gene_list)):
			#TODO: remove the strip when changing the bed reference
			genes_str += "(" + gene_list[i].strip('transcript:') + ":" + str(distance_list[i]) + ")"
		interupeted_gene = "."
		try:
			interupeted_id = distance_list.index(0)
			interupeted_gene = gene_list[interupeted_id]
		except ValueError:
			print("No interupted genes")
		res_row = coord + [interupeted_gene, genes_str]
		res_list.append(res_row)
	res = pd.DataFrame(res_list, columns = res_columns)
	virus = virus.replace(' ', '_')
	res.to_csv(f'{output_directory}/{sample}_{virus}_ins.tsv', sep = '\t', index=False)
	#Output results


#_____________________________________________________________

#load virus list
with open(args.virus_list, 'r') as in_list:
	virus_list = in_list.readlines()
virus_list = [x.strip() for x in virus_list]
print("Virus list")
print(virus_list)

print("\n"+"filtering summary")
filtered_df = filter_summary(args.summary, virus_list, args.o)
print(filtered_df)

print("\n"+"filtering VCF and BAM")

# compress vcf for the downstream tools (tabix)
vcf = args.vcf
cmd = f"bgzip -fc {vcf} > {vcf}.gz"
os.system(cmd)
#Index vcf
vcf = f"{vcf}.gz"
cmd = f"tabix -p vcf {vcf}"
os.system(cmd)

#for each virus detected
if not filtered_df.empty :
	for row in filtered_df.iterrows():
		print("\n"+"filtering VCF")
		print(f'DEBUG:{vcf}')
		print(row[1]['rname'])
		try:
			coord_list = filter_vcf(vcf, row[1]['rname'])
			print(coord_list)
		except ValueError:
			sys.exit()
		print("\n"+"filtering BAM")
		coord_list_updated, sample = fragment_cov_bam(args.bam, coord_list)
		print(coord_list_updated)
		print(sample)
		print("\n"+"getting closest gene")
		closest_genes(coord_list_updated, args.bed, args.o, sample)
		print("\n"+"formating output")
		insertion_report(coord_list_updated, sample, row[1]['name_assigned'], args.o)
		print(insertion_report)
	


# For each line of the virusbreaken summary:

	# Filter res from virus list x

	# Get genome coord and BVF from vcf x

	# Sort, index, filter bam (TO TEST filter bam for both paired match and compare coverage) x

	# Run coverage x

	# Calculate vaf x

	# Get nearby genes (bedtools closest bed) x

	# Output final tables x

	# Link potential start and end of the insertion ?

	# Catch if no insertion are detected

	# decoment bedtool calls