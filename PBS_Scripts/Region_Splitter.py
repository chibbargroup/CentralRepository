'''
Region Splitter, J. Hayes, 3/10/207

Usage: python Region_Splitter.py <fai file> <out_dir>

Purpose: Breaks up a reference genome into several smaller, more manageable regions in order
to speed up the mpileup process. Run mpileup and bcftools on these smaller regions in parallel
in order to enhance process efficiency

Re module used to split up regions based on pattern reconginition (i.e., by chromosome). Change
code starting at line 27 to generate regions as needed. (Change pattern search term to fit needs.)
Will also need to change region list dictionary keys; output file names are the dictionary keys.

Also generates list of region files generated for use with VC_Writer region mode 

'''

from os.path import isfile, isdir, join, basename, dirname
import sys, re

if len(sys.argv) != 3:
	print("whoops, that doesn't work")
else:
	file = sys.argv[1]
	save_dir = sys.argv[2]


region_list = {'AL': [], 'BL': [], 'DL': [], 'AS': [], 'BS': [], 'DS': [], 'U': [], 'B': []}
with open(file, 'r') as f:
	for line in f:
		region_name = line.split('\t')
		if re.match('.*\_\dAL', region_name[0]):
			region_list['AL'].append((region_name[0], region_name[1]))
		elif re.match('.*_\dAS', region_name[0]):
			region_list['AS'].append((region_name[0], region_name[1]))
		elif re.match('.*_\dBL', region_name[0]):
			region_list['BL'].append((region_name[0], region_name[1]))
		elif re.match('.*_\dBS', region_name[0]):
			region_list['BS'].append((region_name[0], region_name[1]))
		elif re.match('.*_\dDL', region_name[0]):
			region_list['DL'].append((region_name[0], region_name[1]))
		elif re.match('.*_\dDS', region_name[0]):
			region_list['DS'].append((region_name[0], region_name[1]))
		elif re.match('.*_\dB', region_name[0]):
			region_list['B'].append((region_name[0], region_name[1]))
		elif re.match('.*_U', region_name[0]):
			region_list['U'].append((region_name[0], region_name[1]))

file_list = []
for key in region_list:
	file_name = join(save_dir, key + '.txt')
	file_list.append(file_name)
	with open(file_name, 'w') as f:
		for region in region_list[key]:
			f.write(region[0] + '\t0\t' + str(int(region[1])-1)  + '\n')

with open(join(save_dir, "region_list.txt"), 'w') as f:
	for file in file_list:
		f.write(file + '\n')