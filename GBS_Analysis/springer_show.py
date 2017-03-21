'''
Springer_Show v2
J. Hayes, Written Mar 21, 2017

Description:
Reads the tab output from the GBS pipeline and assigns each sample allele to either parent A
or parent B for linkage map analysis. The output is a csv file with each allele labelled as 
either "a" (parent a), "b" (parent b), "h" (heterozygous), or "-" (no assignment).

Note: There is no search implementation for column labels in the software, it is therefore necessary
to relabel the columns in the .tab file to remove the file paths and extensions using Excel or a similar
program.

Usage:
python springer_show.py <input file>

User Input: While running, the program will as the user for the labels of Parent A1, Parent A2,
Parent B1, and Parent B2; these must match exactly the column labels in the input .tab file

Output:
Outputs assigned_alleles.csv in the same folder as the as the input data file

'''

import sys, time, argparse
import pandas as pd
from collections import OrderedDict, Counter
from os.path import isfile, isdir, join, dirname
from os import mkdir

##MODIFY THIS TO ADD RE EXPRESSION SEARCH FOR FILE NAME
def File_Reader(file):
	print("Opening file...")
	start_time = time.time()
	df = pd.read_csv(file, delimiter = '\t', header = 0)	
	print("Relabelling...")
	relabelled_index = []
	for row in df.index:
		new_index = "%s_%s" %(df.ix[row, 0], df.ix[row, 1])
		if new_index not in relabelled_index:
			relabelled_index.append(new_index)
		else:
			relabelled_index.append(new_index + '_1')
	df.index = relabelled_index
	df = df.drop(df.columns[0:3], axis = 1)
	print("--- %s seconds ---" % (time.time() - start_time))
	return df

def Parent_Parse(df, parent_dict):
	print("Parsing out the parents...")
	start_time = time.time()
	parent_a = df[parent_dict["A"]]
	parent_b = df[parent_dict["B"]]
	parent_df = pd.concat([parent_a, parent_b], axis = 1)
	print("--- %s seconds ---" % (time.time() - start_time))
	return parent_df

def Parent_Genotyper(df):
	print("Assigning parent genotypes...")
	start_time = time.time()
	genotype_dict = {}
	identical_sites = {}
	non_matching_sites = {"A": {}, "B": {}}
	for location in df.index:
		parent_a1 = df.ix[location, 0]
		parent_a2 = df.ix[location, 1]
		parent_b1 = df.ix[location, 2]
		parent_b2 = df.ix[location, 3]
		if parent_a1 == parent_a2 and parent_b1 == parent_b2 and parent_a1 != './.' and parent_b1 != './.':
			if parent_a1 != parent_b1:
				genotype_dict[location] = {"A": parent_a1, "B": parent_b1}
			else:
				identical_sites[location] = {"A": parent_a1, "B": parent_b1}
		elif parent_a1 != parent_a2:
			non_matching_sites["A"][location] = [parent_a1, parent_a2]
		elif parent_b1 != parent_b2:
			non_matching_sites["B"][location] = [parent_b1, parent_b2]
	print("--- %s seconds ---" % (time.time() - start_time))
	return genotype_dict

def Sample_Genotyper(df, parent_assignments):
	print("Assigning sample genotypes...")
	start_time = time.time()
	df = df.ix[parent_assignments.keys()]
	
	for location in parent_assignments.keys():
		for sample in df.columns:
			if parent_assignments[location]["A"] == df.ix[location, sample]:
				df.set_value(location, sample, "a")
			elif parent_assignments[location]["B"] == df.ix[location, sample]:
				df.set_value(location, sample, "b")
			else:
				assignment = Check_Heterozygous(df.ix[location, sample], parent_assignments[location])
				df.set_value(location, sample, assignment)
	print("--- %s seconds ---" % (time.time() - start_time))
	return df

def Check_Heterozygous(sample_allele, parent_alleles):
	parent_a = parent_alleles["A"]
	parent_b = parent_alleles["B"]

	hetero_1 = parent_a.split('/')[0] + '/' + parent_b.split('/')[0]
	hetero_2 = parent_a.split('/')[1] + '/' + parent_b.split('/')[0] 
	hetero_3 = parent_a.split('/')[0] + '/' + parent_b.split('/')[1]
	hetero_4 = parent_a.split('/')[1] + '/' + parent_b.split('/')[1]
	hetero_list = [hetero_1, hetero_2, hetero_3, hetero_4]

	if sample_allele in hetero_list:
		return "h"
	else:
		return "-"

def Parent_Label_Input():
	correct = False
	yes = ['y', 'yes', 'yup']
	no = ['n', 'no', 'nope']
	while not correct:
		parent_a1 = input("Please input the sample label for Parent A1: ")
		parent_a2 = input("Please input the sample label for Parent A2: ")
		parent_b1 = input("Please input the sample label for Parent B1: ")
		parent_b2 = input("Please input the sample label for Parent B2: ")
		parent_dict = {'A': [parent_a1, parent_a2], 'B': [parent_b1, parent_b2]}
		print("Here are the parents you input: ")
		print(parent_dict)
		correct = Yes_No_Response("Is this correct? [y/n]: ")
	return parent_dict

def Yes_No_Response(question_string):
	positive_response = ['yes', 'y', 'yup', 'yea', 'yeah']
	negative_response = ['no', 'n', 'nope', 'nah', 'nay']
	valid_response = False
	while not valid_response:
		response = input(question_string)
		if response.lower() in positive_response:
			valid_response = True
			return True
		elif response.lower() in negative_response:
			valid_response = True
			return False
		else:
			print("I'm sorry Dave, I can't do that. Please input a valid response.")

def __Main__():
	if len(sys.argv) < 2:
		print("Sorry, you need add the input file as an argument when calling me...")
	elif len(sys.argv) > 2:
		print("Sorry, there are too many arguments, please only give me the input file to analyze")
	else:
		file = sys.argv[1]
		file_dir = dirname(file)
		parent_dict = Parent_Label_Input()
		df = File_Reader(file)
		parent_df = Parent_Parse(df, parent_dict)
		parent_genos = Parent_Genotyper(parent_df)
		final_df = Sample_Genotyper(df, parent_genos)
		final_df.to_csv(join(file_dir, "assigned_alleles.csv"))

__Main__()