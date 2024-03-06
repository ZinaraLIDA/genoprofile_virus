import argparse
import subprocess
import os

def readCoordList(file_path):
	coord_list = []
	try:
		with open(file_path, 'r') as file:
			for line in file:
				if line != "" :
					# Split the tab-separated values into a list
					data = line.strip().split('\t')
					coord_list.append(data)
	except FileNotFoundError:
		print(f"File '{file_path}' not found.")
	except Exception as e:
		print(f"An error occurred: {e}")
	return coord_list

def get_closest_genes(updated_coord_list, bed, output_directory):
	#uses bedtools closest to find genes around the insersion
	#Create bed file for the insertion
	sample = args.updated_coord_list.split("/")[-1].replace('_updated_coord_list.txt', '')
	insertion_bed_file = f'{output_directory}/{sample}.insertion.bed'
	updated_coord_list = readCoordList(updated_coord_list)
	with open(insertion_bed_file, 'w') as out_bed:
		for elem in updated_coord_list:
			print(elem)
			chrom = elem[3]
			pos = elem[4]
			bed_line = [chrom, str(pos), str(int(pos) + 1)]
			out_bed.write('\t'.join(bed_line) + '\n')
			print(bed_line)
	#intersect with bed tools
	output_file = f'{output_directory}/{sample}.gene_distance.tsv'
	cmd = f'bedtools closest -a {insertion_bed_file} -b {bed} -k 5 -D ref > {output_file}'
	os.system(cmd)
	print('bedtools done')

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--updated_coord_list", type=str, help="Path to the updated coordinate list updated")
	parser.add_argument("--bed", type=str, help="Path to the BED file")
	parser.add_argument("--output_directory", type=str, help="Path to the results output directory")
	args = parser.parse_args()
	
	get_closest_genes(args.updated_coord_list, args.bed, args.output_directory)