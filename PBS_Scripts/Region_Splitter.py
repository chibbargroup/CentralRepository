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
		region_name = line.split('\t')
		if re.match('morex', region_name[0]):
			region_list['morex'].append((region_name[0], region_name[1]))
		elif re.match('\d.', region_name[0]):
			region_list['H-chrom'].append((region_name[0], region_name[1]))
		elif re.match('\d', region_name[0]):
			region_list['chrom'].append((region_name[0], region_name[1]))

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