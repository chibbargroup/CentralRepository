from Bio import SeqIO, Entrez
import pandas as pd 
import numpy as np
from os.path import join, isfile, isdir
from os import remove, mkdir, listdir
import sys
Entrez.email = 'john.hayes@usask.ca'

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

#Set the path and list of files to analyze
def GenBankParse():
	path = './GenBank/'
	all_files = [i for i in listdir(path) if isfile(join(path,i))]

	#Put all data into a single dictionary
	data_dict = {}
	for accession in all_files:
		#Open each GB file
		gb_file = join(path, accession)
		record = next(SeqIO.parse(gb_file, 'gb'))
		
		#Pull out some basic info from the record
		key = accession.strip('.gb')
		organism = record.annotations['organism']
		description = record.description
		
		#Parse through the taxonomy list, pick out useful values
		tax = record.annotations['taxonomy']
		taxonomy = ''
		read = False
		for item in tax:
			if item == 'Magnoliophyta':
				read = True
			if read:
				taxonomy += (item + ', ')
		taxonomy = taxonomy.strip(',')

		#Grab the taxonomy id number
		for feature in record.features:
			if feature.type == 'source':
				tax_id = feature.qualifiers['db_xref'][0]
				tax_id = tax_id.strip('taxon:')
		
		data_dict[key] = [organism, tax_id, description, taxonomy]

	#Transfer dictionary into data frame, which then gets saved as CSV file
	df = pd.DataFrame.from_dict(data_dict, orient = 'index')
	df.to_csv('Combined_List.csv', header = False)

def CSV_Parser(file):
	df = pd.read_csv(file)
	accession_list = list(df['Accession Number'])
	holder = []
	for accession in accession_list:
		if accession not in holder:
			holder.append(accession)
	return holder

def File_Reader(file):
	with open(file, 'r') as f:
		i = 1
		holder = []
		for line in f:
			line = line.strip('\n')
			if len(line) > 0:
				if i % 3 == 0:
					accession = line.split(' ')[0]
					if accession not in holder:
						holder.append(accession)
				i += 1
	return holder			

#Reads the results from a blast search, and puts the accession numbers into a single list
def Blast_Parser(file):
	df = pd.read_csv(file, header = 7)
	holder = list(df['Subject_ID'])
	accessions = []

	for item in holder:
		item = item.split(';')
		item = item[0].split('|')
		item = item[3]
		if item not in accessions:
			accessions.append(item)

	return accessions

#Takes two lists, a and b, and combines them such that there are no repeats between the two
#Combines list_b into list_a, and returns list_a with the added results
def List_Combiner(list_a, list_b):
	for item in list_b:
		if item not in list_a:
			list_a.append(item)

	return list_a

#Grabs the genbank information for each entry and saves them as a text file
def NCBI_Genbank_Grabber(accession_list, save_path):
	if not isdir(save_path):
		mkdir(save_path)

	for accession in accession_list:
		file_name = accession + '.gb'
		if not isfile(join(path, file_name)):
			with open(join(path, file_name), 'w') as f:
				print("Fetching %s" %accession)
				handle = Entrez.efetch(db = 'protein', id=accession, rettype = 'gb', retmode = 'text')
				f.write(handle.read())
		else:
			print("The entry %s has already been fetched" %accession)

