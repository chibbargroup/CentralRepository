from os.path import isfile, isdir, join, basename, dirname
import sys, re

if len(sys.argv) != 3:
	print("whoops, that doesn't work")
else:
	file = sys.argv[1]
	save_dir = sys.argv[2]


region_list = {'morex': [], 'H-chrom': [], 'chrom': []}
with open(file, 'r') as f:
	for line in f:
		region_name = line.split('\t')[0]
		if re.match('.?_/dAL', region_name):
			region_list['AL'].append(region_name)
		elif re.match('.?_\dAS', region_name):
			region_list['AS'].append(region_name)
		elif re.match('.?_\dBL', region_name):
			region_list['BL'].append(region_name)
		elif re.match('.?_\dBS', region_name):
			region_list['BS'].append(region_name)
		elif re.match('.?_\dDL', region_name):
			region_list['DL'].append(region_name)
		elif re.match('.?_\dDS', region_name):
			region_list['DS'].append(region_name)
		elif re.match('.?_\dB', region_name):
			region_list['B'].append(region_name)
		elif re.match('_U', region_name):
			region_list['U'].append(region_name)

file_list = []
for key in region_list:
	file_name = join(save_dir, key + '.txt')
	file_list.append(file_name)
	with open(file_name, 'w') as f:
		for region in region_list[key]:
			f.write(region + '\t1\n')

with open(join(save_dir, "region_list.txt"), 'w') as f:
	for file in file_list:
		f.write(file + '\n')