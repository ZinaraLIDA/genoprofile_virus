import argparse
import pandas as pd

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

def insertion_report(updated_coord_list, closest_tsv, output_directory):
	# Define column names
	columns = ['ins_chr', 'ins_start', 'ins_start1', 'an_chr', 'an_start', 'an_end' , 'an_id', 'score', 'strand', 'source', 'feature', 'frame', 'attributes', 'dist']
	res_columns = ['Sample', 'rname', 'virus', 'Chr', 'Pos', 'supporting_frag', 'host_frag', 'vaf', 'interupeted_gene', 'nearby_genes']
	# Read updated_coord_list file
	updated_coord_list = readCoordList(updated_coord_list)
	# Extract sample name from the updated_coord_list path
	sample = args.updated_coord_list.split("/")[-1].replace('_updated_coord_list.txt', '')
	# Load gene distance data from TSV file
	distance_df = pd.read_csv(closest_tsv, sep = '\t', header = None, names = columns)
	print("----------------distance_df------------------")
	print(distance_df)
	print("---------------dimension distance_df--------------------")
	print(distance_df.shape)
	print("---------------head distance_df-------------------------")
	print(distance_df.head())
	for coord in updated_coord_list:
		print("---------------coord-------------------------")
		print(coord)
		sub_df = distance_df[(distance_df["ins_chr"].astype('str') == str(coord[3])) &
		 (distance_df["ins_start"].astype('str') == str(coord[4]))]
		print("-----------------sub_df------------------")
		print(sub_df)
		gene_list = sub_df['an_id'].to_list()
		distance_list = sub_df['dist'].to_list()
		genes_str =''
		for gene, distance in zip(gene_list, distance_list):
			genes_str += f"({gene.strip('transcript:')}:{str(distance)})"
		interupeted_gene = "."
		try:
			interupeted_id = distance_list.index(0)
			interupeted_gene = gene_list[interupeted_id]
		except ValueError:
			print("No interupted genes")
		res_row = coord + [interupeted_gene, genes_str]
		print("---------------------res_row-------------------")
		print(res_row)
		
		# Create a Pandas DataFrame for the result
		res = pd.DataFrame(columns = res_columns)
		res.loc[len(res)] = res_row
		# Replace spaces with underscores in the virus name
		virus = res_row[2]
		virus = virus.replace(' ', '_')
		# Save the insertion report to a TSV file
		res.to_csv(f'{output_directory}/{sample}_{virus}_ins.tsv', sep = '\t', index=False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Generate insertion report.")
	parser.add_argument("--updated_coord_list", type=str, required=True, help="Path to the filtered coordinate list")
	parser.add_argument("--closest_tsv", type=str, required=True, help="Path to the gene_distance file")
	parser.add_argument("--output_directory", type=str, required=True, help="Output directory")
	args = parser.parse_args()

	insertion_report(args.updated_coord_list, args.closest_tsv, args.output_directory)