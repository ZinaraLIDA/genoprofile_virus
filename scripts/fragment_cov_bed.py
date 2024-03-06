import subprocess
import pandas as pd
import argparse
import re

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

def run_bedtools(input_file, bedtools_output_file):
	cmd = f'bedtools genomecov -ibam {input_file} -pc -bga > {bedtools_output_file}'
	subprocess.call(cmd, shell=True)

def fragment_cov_bed(coord_list, samtools_output_file, output_directory, window_size=4000):
	sample_name = samtools_output_file.split("/")[-1]
	print(f"sample_name : {sample_name}")
	#pattern = r"(.+?)_\d+_\d+_\d+\.bam"
	pattern = re.compile(r'([0-9A-Z]+_[A-Z-]+)')
	#match = re.search(pattern, sample_name)
	match = pattern.search(sample_name)
	print(f"match : {match}")
	sample_name = match.group(1)
	updated_coord_list = []
	coord_list = readCoordList(coord_list)
	for coord in coord_list:
		rname, virus, chrom, pos, bvf = coord
		start = int(int(pos) - (window_size // 2))
		end = int(int(pos) + (window_size // 2))
		bedtools_output_file = f'{output_directory}/{sample_name}_{chrom}_{start}_{end}.bed'
		run_bedtools(samtools_output_file, bedtools_output_file)
		cov_df = pd.read_csv(bedtools_output_file, sep='\t', header=None, names=['chr', 'start', 'end', 'cov'])
		# Drop the last row which is 'None'
		cov_df = cov_df.dropna()
		cov_df = cov_df.astype({'start': int, 'end': int, 'cov': int})
		# Find the row that matches the position and chromosome
		cov_row = cov_df[(cov_df['chr'] == chrom) & (cov_df['start'] <= int(pos)) & (cov_df['end'] > int(pos))]
		if not cov_row.empty:
			cov = int(cov_row.iloc[0]['cov'])
			updated_coord = [rname, virus, chrom, pos, bvf, cov, int(bvf) / (int(bvf) + cov)]
		else:
			# Handle the case where the position is not found in the coverage data
			updated_coord = [rname, virus, chrom, pos, bvf, 0, 1.0]
		updated_coord_list.append(updated_coord)
	# Remove 'X' and 'Y' chromosomes from the list
	updated_coord_list = [coord for coord in updated_coord_list if coord[2] not in ['X', 'Y']]
	updated_coord_list.sort()
	sample_name = sample_name.strip(".bam")
	return updated_coord_list, sample_name

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--coord_list", type=str, help="Path to the coord_list file")
	parser.add_argument("--samtools_output_file", type=str, help="Path to the samtools output file")
	parser.add_argument("--output_directory", type=str, help="Directory for output")
	args = parser.parse_args()
	
	updated_coord_list, sample_name = fragment_cov_bed(args.coord_list, args.samtools_output_file, args.output_directory)

	sample = args.coord_list.split("/")[-1].replace('_coord_list.txt', '')
	with open(f"{args.output_directory}/{sample}_updated_coord_list.txt", "w") as fileout :
		for coord in updated_coord_list:
			print("----------------coord-------------")
			print(coord)
			fileout.write(f"{sample_name}"+"\t"+f"{coord[0]}"+"\t"+f"{coord[1]}"+"\t"+f"{coord[2]}"+"\t"+f"{coord[3]}"+"\t"+f"{coord[4]}"+"\t"+f"{coord[5]}"+"\t"+f"{coord[6]}"+"\n")
