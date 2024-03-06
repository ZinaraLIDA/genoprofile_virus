import argparse
import pandas as pd

def readVirusList(virusListPath):
    with open(virusListPath, 'r') as in_list:
        virus_list = in_list.readlines()
        virus_list = [x.strip() for x in virus_list]
    print(virus_list)
    return virus_list

def filter_summary(summary, virus_list, output_directory):
    summary_df = pd.read_csv(summary, sep='\t', header=0)
    sample = summary.split("/")[-1].replace('.vcf.summary.tsv', '')
    virus_list = readVirusList(virus_list)
    filtered_df = summary_df[(summary_df["name_species"].isin(virus_list)) | (summary_df["name_genus"].isin(virus_list))]
    if not filtered_df.empty:
        filtered_df.to_csv(f"{output_directory}/{sample}_all_virus.tsv", sep='\t', index=False)
    else:
        filtered_df = pd.DataFrame()
        with open(f"{output_directory}/{sample}_all_virus.tsv", "w") as filin:
            pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", type=str, help="Path to the summary file")
    parser.add_argument("--virus_list", type=str, help="Path to the virus list")
    parser.add_argument("--output_directory", type=str, help="Output directory")
    args = parser.parse_args()
    
    filter_summary(args.summary, args.virus_list, args.output_directory)