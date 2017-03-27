from os.path import isdir, dirname
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
		' for cases where there are a large number of SNPs; Default is False')
	parser.add_argument('--exclude_heterozygotes', '-h', type = bool, default = False, help = 'Choose if you want to exclude sites with' + \
		' heterozygous assignments (h); Default is True')
	return parser

def __Main__():
	parser = ArgParse_Init().parse_args()
	if not isdir(dirname(parser.out_file)):
		print("Creating the directory %s" %dirname(parser.outfile))
		mkdir(dirname(parser.out_file))
	if parser.remove_duplicates:
		Loci_Relabeller(parser.in_file, parser.out_file)
		Scaffold_Repeat_Remover(parser.out_file, parser.out_file)
	if parser.output_type == 'mst':
		
	elif parser.output_type == 'joinmap':
		print("joinmap selected")

def Read_Tab_File(in_file):
	file_contents = OrderedDict()
	with open(in_file, 'r') as f:
		f.readline()
		for line in f:
			split_line = line.split('\t')
			locus_label = split_line[0]
			data = '\t'.join(split_line[1:])
			file_contents[locus_label] = data
	return file_contents

def Loci_Relabeller(data):
	relabeled_data = OrderedDict()
	for key in data:
		split_locus = key.split('_')
		new_locus = split_locus[2] + '_' + split_locus[3]
		i = 1		
		#Deal with case where there are multiple loci with same name
		while new_locus in new_data.keys():
			new_locus = new_locus + '_%s' %i
			i += 1
		relabelled_data[new_locus] = data[key]
	return relabelled_data

def Scaffold_Repeat_Remover(data, remove_h):
	scaffolds = OrderedDict()
	with open(in_file, 'r') as f:
		header = f.readline()
		for line in f:
			split_line = line.split('\t')
			blank_count = split_line.count('-')
			scaffold_number = split_line[0].split('_')[2] + split_line[0].split('_')[3]
			if 'h' in split_line[1:-1]:
				pass
			elif scaffold_number not in scaffolds.keys():
				scaffolds[scaffold_number] = (blank_count, line)
			else:
				if blank_count < scaffolds[scaffold_number][0]:
					scaffolds[scaffold_number] = (blank_count, line)
	return scaffolds

def JM_Header(parameters):
	header = "name=%s\n" %parameters[0] + \
	"popt=%s\n" %parameters[1] + \
	"nloc=%s\n" %parameters[2] + \
	"nind=%s\n\n" %parameters[3]
	return header

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
	"number_of_individual %s\n" %parameters[11]
	return header

def File_Writer(scaffolds, out_file):	
	with open(out_file, 'w') as f:
		f.write(header)
		for key in scaffolds:
			f.write(scaffolds[key][1])

def MST_Parameter_Input():
	pass

def JM_Paramter_Input():
	pass

Main()


