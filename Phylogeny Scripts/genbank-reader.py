'''
GenBank Reader
Written by J. Hayes; last edit 12/9/2016

Purpose: Takes a csv file of protein or DNA sequence accession numbers and names and reads the corresponding
genbank file (assuming it's already been grabbed). Then searches for features in the GenBank file and outputs
these features in a fasta file for further analysis

CSV File should have the following column headers: Accession Number, Common Names, Familiy, Protein Name (case sensitive)
Usage:
python genbank-reader.py csv_file_loc feature_type feature_id genbank_file_dir output_dir

feature_type: The type of feature to be retrieved (i.e., region, CDS, site, etc)
feature_id: Some sort of string that will help the program the feature name and pick out the correct one

Output:
Saves a file entitled sites.fa in the output directory

NOTE: Will rewrite file if sites.fa file already exists
'''

from Bio import SeqIO
import pandas as pd 
import numpy as np
from os.path import join, isfile
from os import remove
import sys

#Reads a CSV file containing accession number, species name, family group, and protein type info
#and returns a dictionary with the accession number as the key; for now, family is based on the 
def CSV_Parser(csv_file):
	df = pd.read_csv(csv_file)
	data_dict = {}
	for row in df.iterrows():
		key = row[1]['Accession Number']
		if pd.notnull(key):
			species = row[1]['Common Names']
			family = row[1]['Family']
			protein_name = row[1]['Protein Name']
			data_dict[key] = [species, family, protein_name]
	return data_dict

#Takes the data dictionary, the sequence feature type (i.e., CDS, Region, Site, etc), some string to identify the region
#by the note (i.e., "ADP Binding Site"), and the file path where the GB files are stored; returns data dictionary with the 
#feature names and locations appended as tuples
def Genbank_Parser(data_dict, seq_type, region_ident, gb_path):
	for key in data_dict:
		#Open the genbank file
		gb_file = join(gb_path, key + '.gb')
		record = next(SeqIO.parse(gb_file, 'gb'))
		
		#Read the genbank file and identify relevant features/regions
		full_sequence = record.seq
		feat_holder = []
		for feature in record.features:
			try:
				if feature.type == seq_type and region_ident in feature.qualifiers['note'][0]:
					feature_name = feature.qualifiers['note'][0]
					feature_loc = feature.location
					holder = []
					for position in feature_loc:
						holder.append(full_sequence[position])
					tup = (feature_name, holder)
					feat_holder.append(tup)
			except KeyError:
				print("One or more regions types of interest were found, but did not have a note attached")
		if len(feat_holder) > 0:
			data_dict[key].append(feat_holder)
		else:
			data_dict[key].append('EMPTY')
	return(data_dict)

#Takes the info from the data dictonary and puts it into a fasta file in the directory specified by the user
#with a separate entry per site
def Site_Writer(data_dict, out_path):
	#Determine if site.fa already exsists; if it does, remove it and let the user it was removed
	if isfile(join(out_path, 'sites.fa')):
		remove(join(out_path, 'sites.fa'))
		print("Caution, I detected a previous site.fa file in this location and deleted it")

	#Write the site information to the file
	with open(join(out_path, 'sites.fa'), 'a') as f:
		for key in data_dict:
			if data_dict[key][-1] != "EMPTY":
				for feature,sequence in data_dict[key][-1]:
					#Write fasta header; format: >species, family, protein name, feature name
					species = data_dict[key][0].strip('\n')
					family = data_dict[key][1].strip('\n')
					protein_name = data_dict[key][2].strip('\n')
					feature = feature.strip('\n')
					f.write('>' + species + ' ' + family + ' ' + protein_name + ' ' + feature + '\n')
					
					#Write the sequence to the file; 70 characters per line
					i = 1
					for char in sequence:
						if i % 70 == 0:
							f.write('\n' + char)
							i += 1
						else:
							f.write(char)
							i += 1
					f.write('\n\n')

#General functions to call and run the above functions in one easy step
def Program_Runner(csv_file, feature_type, feature_id, gb_file_path, out_path):
	data_dict = CSV_Parser(csv_file)
	data_dict = Genbank_Parser(data_dict, feature_type, feature_id, gb_file_path)
	Site_Writer(data_dict, out_path)

'''
#Temporary parameters for bug testing
csv_file = './Protein list-Branching.csv'
feat_type = 'Region'
feature_id = 'glucan'
gb_file_path = './Protein_sequence/Branching/Genbank'
out_path = './'
Program_Runner(csv_file, feat_type, feature_id, gb_file_path, out_path)
'''

#Check the lenght of sys.argv; it the appropriate number of arguments are available, run the program!
if len(sys.argv) < 6:
	print('Whoops, not enough arguments; current number of args is %i' %len(sys.argv))
	print('I need the following (in order): List location, Feature Type, Feature ID, Genbank file path, and output directory')
elif len(sys.argv) > 6:
	print('Whoops, not enough arguments; current number of args is %i' %len(sys.argv))
	print('I need the following (in order): List location, Feature Type, Feature ID, Genbank file path, and output directory')
else:
	Program_Runner(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])



