from os.path import join
import sys

region_file = ".\chome_list.txt"
ref_dir = './References'
save_dir = './output'

def MPileup_File_Template(ref_location, out_name):
	file_template = """###
# The PBS Header Portion (Resource Request Settings)
#
###

##Set the shell
#! /bin/bash
#PBS	-S	/bin/bash		

##Main request settings: no. of processors, walltime, ram request
#PBS	-l	procs=1		
#PBS	-l	walltime=72:00:00
#PBS	-l	pmem=10gb	

##Email settings; a = when aborted, b = when started, e = when finished
#PBS	–M	d.sanche14@gmail.com
#PBS	-m	abe		)

###
#
# The Script Portion
#
###
## Initialize Programs
module load application/samtools/1.3.1
cd $PBS_O_WORKDIR
cd ../data

## Run Program
samtools mpileup -Duf %s -b ./bam_list.txt > %s.raw.bcf""" %(ref_location, out_name)
	return file_template

region_list = []
with open(region_file, 'r') as f:
	for line in f:
		line = line.strip('\n')
		region_list.append(line)

for region in region_list:
	ref_location = join(ref_dir, region)
	region_info = region.split('.')
	region_name = region_info[-3] + "_" + region_info[-2]
	pbs_file_name = join(save_dir, region_name + '.pbs')
	file = MPileup_File_Template(join(ref_dir, region), region_name)
	with open(pbs_file_name, 'w') as f:
		f.write(file)
