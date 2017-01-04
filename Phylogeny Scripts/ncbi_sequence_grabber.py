'''
Script "suite" for grabbing and analyzing files from the NCBI database. Generally, the program will
take a list of accession numbers and retrieve either the protein or nucleotide sequence and can save the
file in either genbank or fasta form. The file containing the accession numbers should be in csv format
with a column labelled as "Accession". Other column labels are also possible, and up to the user.

Usage:
python ncbi_sequence_analyzer.py 'accession list', 'type to retrieve'*, 'path to save files'
*type to retrieve can equal 'protein' or 'nucleotide'

'''
from Bio import Entrez
from os.path import join, basename, isfile
from os import listdir
import pandas as pd
import numpy as np 
import sys

#Always tell NCBI your email
Entrez.email = 'john.hayes@usask.ca'

#Fetches either individual protein or nucleotide entry from database and saves it in specified path
#Default format fetched in genbank; entry_type = 'protein' or 'nucleotide'
def Sequence_Fetch(entry_id, entry_type, format_type, path):
	#Format file name appropriately
	if format_type == 'gb':
		file_name = entry_id + '.gb'
	elif format_type == 'fasta':
		file_name = entry_id + '.fasta' 
	
	#Fetch the sequence and sav eit
	with open(join(path, file_name), 'w') as f:
		print("Fetching %s" %entry_id)
		handle = Entrez.efetch(db = entry_type, id=entry_id, rettype = format_type, retmode = 'text')
		f.write(handle.read())

#Reads the list of accession numbers. The file containing the accession numbers should be in csv format
#with a column labelled as "Accession Number"
def Accession_Reader(accession_file, entry_type, format_type, path):
	#Set the extension based on the file type
	if format_type == 'fasta':
		extension = '.fasta'
	elif format_type == 'genbank' or format_type == 'gb':
		extension = '.gb'
	else:
		print("Error, not a valid file type")
		return

	#Not necessary to read as dataframe, but likely easier to do this for future file name/fasta header manipulation
	#Use of SeqIO and parsing through genbank files may make this unnecessary
	df = pd.read_csv(accession_file)
	
	#Make a list of the accession numbers to feed to the sequence fetcher
	accession_list = []
	for entry in list(df['Accession Number']):
		if pd.notnull(entry):
			accession_list.append(entry)
	
	for entry in accession_list:
		if not isfile(join(path, entry + extension)):
			Sequence_Fetch(entry, entry_type, format_type, path)
		else:
			print("The file with Accession Number %s appears to have already been fetched" %entry)

#Run the program
if len(sys.argv) == 4:
	Accession_Reader(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
else:
	print("Sorry, I don't follow, there are supposed to be four arguments...")
	print("'accession list', 'type to retrieve', format type, 'path to save files'")




