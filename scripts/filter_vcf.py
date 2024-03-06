import argparse
import re
import pandas as pd
from pysam import VariantFile

def filter_vcf(vcf, filtered_summary):
	coord_list = []
	in_vcf = VariantFile(vcf)
	try :
		filtered_summary = pd.read_csv(filtered_summary, sep = '\t', header = 0)
		for row in filtered_summary.iterrows():
			print("----------------row-------------")
			print(row)
			rname = row[1]['rname']
			print("----------------rname-------------")
			print(rname)
			virus = row[1]['name_assigned']
			print("----------------virus-------------")
			print(virus)
			for rec in in_vcf.fetch(rname):
				print("--------------------rec-----------------")
				print(rec)
				bvf = rec.info["BVF"]
				print(f"bvf : {bvf}")
				alt = rec.alts[0]
				print(f"alt : {alt}")
				pos_re = r"([0-9]+|[X,Y]):([0-9]+)"
				m = re.search(pos_re, alt)
				print(f"m : {m}")
				chrom = m.group(1)
				print(f"chrom : {chrom}")
				pos = m.group(2)
				print(f"pos : {pos}")
				coord_list.append([rname, virus, chrom, pos, bvf])
	except pd.errors.EmptyDataError:
		print('CSV file is empty')
	except FileNotFoundError:
		print('CSV file not found')
	print("--------------coord_list------------")
	print(coord_list)
	return coord_list

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--vcf", type=str, help="Path to the VCF file")
	parser.add_argument("--filtered_summary", type=str, help="Summary file filtered")
	parser.add_argument("--output_directory", type=str, help="Directory for output")
	args = parser.parse_args()
	
	updated_coord_list = filter_vcf(args.vcf, args.filtered_summary)
	
	sample = args.filtered_summary.split("/")[-1].replace('_markdup_all_virus.tsv', '')
	with open(f"{args.output_directory}/{sample}_coord_list.txt", "w") as fileout :
		for coord in updated_coord_list:
			print("----------------coord-------------")
			print(coord)
			fileout.write(f"{coord[0]}"+"\t"+f"{coord[1]}"+"\t"+f"{coord[2]}"+"\t"+f"{coord[3]}"+"\t"+f"{coord[4]}"+"\n")