'''
Genbanker.py
Written by J. Hayes, last edit 1/10/17

Purpose: Reads through lists of accession numbers generated through similarity comparisons,
blast searches, and general NCBI searches and fetches the corresponding Genbank file

Input: Blast search (.csv; file must have 'blast' or 'Blast' in the file name)
	   Accession list exported from web database search
	   Similarity files (.csv format)
	   Any csv file with accession numbers in a column; column should be labelled as 'axcn' or
	   change line ~67 in General_CSV_Parser function as needed

Output: Saves the genbank files in the specified output directory as .gb files; format of name is
		accession_number.gb

Usage: python Genbanker.py Input_File_Directory Save_Directory Retrieval_Type
	*Note: Retrieval type can be any of the following: Protein, protein, prot, Nucleotide, nucleotide, nt
	(Other options to specifiy protein or nucleotide can be added to the lists in lines ~127, 128 in Program Runner)

'''

from os.path import isfile, isdir, join
from os import mkdir, listdir
from Bio import Entrez
import pandas as pd
import sys

### Always tell NCBI your email; change when used by a different user
Entrez.email = 'john.hayes@usask.ca'

#Parses through the blast results table and pulls out the GB accession numbers; function is based 
#on the output from the NCBI Blast results table export; returns the accession numbers in a list
def Blast_Parser(file):
	#Read the CSV; generally the data headers are rows 1-7, and the data columns are in row 8
	#(Row 8 indexes to row 7 bc Python starts at 0)
	df = pd.read_csv(file, header = 7)
	subject_ids = list(df['Subject_ID'])
	holder = []
	gb_accessions = []

	#Split out the results; in some cases there are multiple accession numbers per row; these entries
	#are separated by a ';'
	for item in subject_ids:
		result = item.split(';')
		for entry in result:
			holder.append(entry)
	
	#Once match sets are split, read the match sets; pick up when a code is 'gb', and capture those
	#codes (parse through and ignore the GI indicies)
	for item in holder:
		result = item.split('|')	
		capture = False
		for entry in result:
			if entry == 'gb':
				capture = True

			if capture and entry != 'gb':
				gb_accessions.append(entry)
				capture = False
	return gb_accessions

#Read the accession numbers from a CSV file with a column labelled "Accession Number"; good
#for results from similarity searches or more general lists
def General_CSV_Parser(file):
	df = pd.read_csv(file)
	accession_list = list(df['acxn']) #Note, the column name may have to be changed based on the input file
	holder = []
	for accession in accession_list:
		if accession not in holder:
			holder.append(accession)
	return holder

#Reads the accession numbers generated from an NCBI search; these search results are exported from the
#NCBI website by clicking Send:; then selecting file, and selecting format: Accession Number, which generates a .txt file
def NCBI_Search_Reader(file):
	accession_list = []
	with open(file, 'r') as f:
		for line in f:
			line = line.strip('\n')
			accession_list.append(line)
	return accession_list		

#Takes two lists, a and b, and combines them such that there are no repeats between the two
#Combines list_b into list_a, and returns list_a with the added results
def List_Combiner(list_a, list_b):
	for item in list_b:
		if item not in list_a:
			list_a.append(item)
	return list_a

#Fetches the protein genbank file for each accession in the list passed to it (accession_list), and saves 
#the genbank files as separate .gb files in the specified file file path (save_path)
def NCBI_Protein_Genbank_Grabber(accession_list, save_path):
	for accession in accession_list:
		file_name = accession + '.gb'
		if not isfile(join(save_path, file_name)):
			with open(join(save_path, file_name), 'w') as f:
				print("Fetching %s" %accession)
				handle = Entrez.efetch(db = 'protein', id=accession, rettype = 'gb', retmode = 'text')
				f.write(handle.read())
		else:
			print("The entry %s has already been fetched" %accession)

#Fetches the nucleotide genbank file for each accession in the list passed to it (accession_list), and saves 
#the genbank files as separate .gb files in the specified file file path (save_path)
def NCBI_Nucleotide_Genbank_Grabber(accession_list, save_path):
	for accession in accession_list:
		file_name = accession + '.gb'
		if not isfile(join(save_path, file_name)):
			with open(join(save_path, file_name), 'w') as f:
				print("Fetching %s" %accession)
				handle = Entrez.efetch(db = 'nucleotide', id=accession, rettype = 'gb', retmode = 'text')
				f.write(handle.read())
		else:
			print("The entry %s has already been fetched" %accession)

def Program_Runner(file_directory, save_directory, rettype):
	#Generate a list of files in the file directory; make the save directory if necessary, and initialize the 
	#accession list 
	file_list = [f for f in listdir(file_directory) if isfile(join(file_directory, f))]
	if not isdir(save_directory):
		mkdir(save_directory)
	accession_list = []
	
	#Lists of options for interpreting rettype; can be added to
	protein = ['Protein', 'protein', 'prot', 'Prot']
	nucleotide = ['Nucleotide', 'nucleotide', 'nt', 'NT', 'Nt']

	#Read each file in the file list; procedure is based on extension and if blast is in the file name
	for file in file_list:
		file = join(file_directory, file)
		if file[-4:] == '.csv' and 'blast' not in file and "Blast" not in file:
			results = General_CSV_Parser(file)
			accession_list = List_Combiner(accession_list, results)
		elif file [-4:] == '.csv':
			results = Blast_Parser(file)
			accession_list = List_Combiner(accession_list, results)
		elif file [-4:] == '.txt':
			results = NCBI_Search_Reader(file)
			accession_list = List_Combiner(accession_list, results)

	#Grab the genbank file for each accession number; use Protein or Nucleotide based on user input when calling the script
	if rettype in protein:
		NCBI_Protein_Genbank_Grabber(accession_list, save_directory)
	elif rettype nucleotide:
		NCBI_Nucleotide_Genbank_Grabber(accession_list, save_directory)
	else:
		print("Retrival type not understood; use Protein or Nucleotide")

if len(sys.argv) > 4 or len(sys.argv) < 4:
	print("I'm sorry, I don't understand; I only take 3 arguments")
	print("These arguments are: 1) File directory, 2) Save directory, 3) Retrieval type")
else:
	if not isdir(sys.argv[1]):
		print("The file directory you gave me doesn't seem to exist...try calling me again with the correct directory")
	else:
		Program_Runner(sys.argv[1], sys.argv[2], sys.argv[3])
