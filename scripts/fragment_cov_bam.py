import subprocess
import pandas as pd
import argparse

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

def run_samtools(bam, chrom, start, end, samtools_output_file):
	cmd = f'samtools view -hb {bam} "{chrom}:{start}-{end}" > {samtools_output_file}'
	subprocess.call(cmd, shell=True)

def fragment_cov_bam(bam, coord_list, output_directory, window_size=4000):
	sample_name = bam.split("/")[-1].strip(".bam")
	updated_coord_list = []
	coord_list = readCoordList(coord_list)
	for coord in coord_list:
		rname, virus, chrom, pos, bvf = coord
		start = int(int(pos) - (window_size // 2))
		end = int(int(pos) + (window_size // 2))
		samtools_output_file = f'{output_directory}/{sample_name}_{chrom}_{start}_{end}.bam'
		run_samtools(bam, chrom, start, end, samtools_output_file)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--bam", type=str, help="Path to the bam file")
	parser.add_argument("--coord_list", type=str, help="Path to the coord_list file")
	parser.add_argument("--output_directory", type=str, help="Directory for output")
	args = parser.parse_args()
	
	fragment_cov_bam(args.bam, args.coord_list, args.output_directory)