from os.path import isdir, dirname, isfile
from os import mkdir
from collections import OrderedDict
import argparse, sys

def ArgParse_Init():
	parser = argparse.ArgumentParser(description="Scrub tab files and convert to .loc files")
	parser.add_argument('--in_file', '-i', type=str, required = True, help = "The input file to be processed")
	parser.add_argument('--out_file', '-o', type=str, required = True, help = "Where to save the output file")
	parser.add_argument('--output_type', '-t', type = str.lower, required = True,
		choices = ['mst', 'joinmap'], help = "MST or joinmap loc file output")
	parser.add_argument('--relabel_loci', '-l', type = bool, required = False, default = False, 
		help = "Use if you need to relabel loci points (currently implemented for wheat sequences only)")
	parser.add_argument('--remove_duplicates', '-r', type = bool, default = False, help = 'Use if you want only one SNP per scaffold' + \
		' for cases where there are a large number of SNPs and scaffolds; Default is False')
	parser.add_argument('--exclude_heterozygotes', '-e', type = bool, default = False, help = 'Choose if you want to exclude sites with' + \
		' heterozygous assignments (h); Default is True')
	return parser

def __Main__():
	parser = ArgParse_Init().parse_args()
	
	#Read the data and make output direcotry if necessary
	if not isfile(parser.in_file):
		print("I can't find that file")
		sys.exit()
	else:
		print("Reading the tab file...")
		data = Read_Tab_File(parser.in_file)
	if not isdir(dirname(parser.out_file)) and dirname(parser.out_file) != "":
		print("Creating the directory %s" %dirname(parser.out_file))
		mkdir(dirname(parser.out_file))
	
	#The Optional Operations
	if parser.exclude_heterozygotes:
		print("Removing heterozygous loci...")
		data = Remove_Heterozygotes(data)
	if parser.remove_duplicates:
		print("Picking a single SNP from each scaffold...")
		data = Scaffold_Repeat_Remover(data)
	if parser.relabel_loci:
		print("Relabelling the loci...")
		data = Loci_Relabeller(data)

	#Required Operations
	if parser.output_type == 'mst':
		print("Looks like we're making an MST Loc File...")
		print("Let me get some info from you...")
		header_params = MST_Parameter_Input(data)
		header = MTS_Header(header_params)
		File_Writer(header, data, parser.out_file)
	elif parser.output_type == 'joinmap':
		print("Looks like we're making a JM Loc File...")
		print("Let me get some info from you...")
		header_params = JM_Parameter_Input(data)
		header = JM_Header(header_params)
		File_Writer(header, data, parser.out_file)

def Read_Tab_File(in_file):
	file_contents = OrderedDict()
	with open(in_file, 'r') as f:
		f.readline()
		for line in f:
			split_line = line.split('\t')
			locus_label = split_line[0] + "_0"
			data = '\t'.join(split_line[1:])
			#Deal with case where there are multiple loci with same name
			i = 1
			while locus_label in file_contents.keys():
				locus_label = locus_label[:-2] + '_%s' %i
				i += 1
			file_contents[locus_label] = data
	return file_contents

def Loci_Relabeller(data):
	relabelled_data = OrderedDict()
	for key in data:
		split_locus = key.split('_')
		new_locus = split_locus[2] + '_' + split_locus[3] + "_0"		
		#Deal with case where there are multiple loci with same name
		i = 1	
		while new_locus in relabelled_data.keys():
			new_locus = new_locus[:-2] + '_%s' %i
			i += 1
		relabelled_data[new_locus] = data[key]
	return relabelled_data

def Remove_Heterozygotes(data):
	new_data = OrderedDict()
	for key in data:
		if 'h' not in data[key].split('\t'):
			new_data[key] = data[key]
	return new_data

def Scaffold_Repeat_Remover(data):
	new_data = OrderedDict()
	blank_count_dict = {}
	for key in data:
		scaffold = '_'.join(key.split('_')[0:-2])
		allele_list = data[key].split('\t')
		sample_blanks = allele_list.count('-')
		if scaffold not in new_data.keys():
			new_data[scaffold] = data[key]
			blank_count_dict[scaffold] = sample_blanks
		else:
			if sample_blanks < blank_count_dict[scaffold]:
				new_data[scaffold] = data[key]
				blank_count_dict[scaffold] = sample_blanks
	return new_data

def JM_Parameter_Input(data):
	name = input("Please input a name for this population: ")
	popt = input("Please input the population type (F2, RIL, DH1, etc.): ")
	nloc = len(data.keys())
	nind = len(list(data.values())[1].split('\t'))
	parameters = [name, popt, nloc, nind]
	return parameters

def JM_Header(parameters):
	header = "name=%s\n" %parameters[0] + \
	"popt=%s\n" %parameters[1] + \
	"nloc=%s\n" %parameters[2] + \
	"nind=%s\n\n" %parameters[3]
	return header

def MST_Parameter_Input(data):
	popt = input("Please input the population type (F2, RIL, DH1, etc.): ")
	pop_name = input("Please input a name for this population: ")
	distance_function = input("Please input the distance function (kosambi or haldane): ").lower
	cut_off_p_value = input("Please input the threshold for clustering (suggested 0.000001): ")
	no_map_dist = input("Please input the minimum mapping distance in centimorgans (suggested 2): ")
	no_map_size = input("Please put the maximum marker separation distance (recommended 15): ")
	missing_threshold = input("Please put the maximum amount of missing markers (as decimal, recommended 0.25): ")
	estimation_before_clustering = input("Estimate missing data before clustering? (yes or no): ").lower
	objective_function = input("What objective function should be used (COUNT or ML)? ").upper
	nloc = len(data.keys())
	nind = len(list(data.values())[1].split('\t'))
	parameters = [popt, pop_name, distance_function, cut_off_p_value, no_map_dist, no_map_size, missing_threshold,
	estimation_before_clustering, objective_function, nloc, nind]
	return parameters

def MTS_Header(parameters):
	header = "population_type %s\n" %parameters[0] + \
	"population_name %s\n" %parameters[1] + \
	"distance_function %s\n" %parameters[2] + \
	"cut_off_p_value %s\n" %parameters[3] + \
	"no_map_dist %s\n" %parameters[4] + \
	"no_map_size %s\n" %parameters[5] + \
	"missing_threshold %s\n" %parameters[6] + \
	"estimation_before_clustering %s\n" %parameters[7] + \
	"detect_bad_data %s\n" %parameters[8] + \
	"objective_function %s\n" %parameters[9] + \
	"number_of_loci %s\n" %parameters[10] + \
	"number_of_individual %s\n\n" %parameters[11]
	return header

def File_Writer(header, data, out_file):	
	with open(out_file, 'w') as f:
		f.write(header)
		for key in data:
			f.write(key + '\t')
			f.write(data[key])

__Main__()


