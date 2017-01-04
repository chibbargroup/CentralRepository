'''
Combines fasta files based on a list of accession numbers with formated headers

Usage:
python sequence_combiner.py file_info.csv fasta_directory out_directory

'''

from os import listdir, remove
from os.path import join, isfile
import pandas as pd 
import sys

def File_Reader(file_info):
	data_dict = {}
	df = pd.read_csv(file_info)

	for row in df.iterrows():
		key = row[1]["Accession Number"]
		name = row[1]["Common Names"]
		family = row[1]["Family"]
		protein = row[1]["Protein Name"]
		data_dict[key] = [name, family, protein]

	return data_dict

def File_Combiner(data_dict, fasta_path, out_path):
	out_file = join(out_path, 'combined.fasta')
	if isfile(out_file):
		remove(out_file)
		print("I've removed the previous file found here, and saved the new file from this run")


	with open(out_file, 'a') as f1:
		for key in data_dict:
			file_name = join(fasta_path, key + '.fasta')
			with open(file_name, 'r') as f2:
				for line in f2:
					if ">" in line:
						f1.write('>' + data_dict[key][0] + ' ' + data_dict[key][1] + ' ' + data_dict[key][2] + '\n')
					else:
						f1.write(line)
			f1.write('\n')

if len(sys.argv) == 4:
	data_dict = File_Reader(sys.argv[1])
	File_Combiner(data_dict, sys.argv[2], sys.argv[3])
else:
	print("Argument length doesn't match; should have 3 arguments:")
	print("1) File.csv, 2) Fasta directory, 3) Output directory")
