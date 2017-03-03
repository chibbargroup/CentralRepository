'''
BT_ST PBS File Writer
Written by J. Hayes, 2-27-2017

Purpose: This script is used to write PBS files to run the bowtie and samtools cleanup pipeline.
It is meant for use with Westgrid, and currently works with Samtools 1.3.1. It will create a PBS file
for .fq file in the specified data directory

NOTE: All file paths must lead to the target from the directory where the PBS files are submitted to Torque

Usage: python BT_ST_writer.py data_dir pbs_dir results_dir ref_genome email (optional)

List of inputs
1) data_dir: Directory containing the data to be aligned
2) pbs_dir: Directory where the created PBS files will be saved
3) results_dir: Directory where analysis results should be saved
4) ref_genome: Path to the bowtie indexed reference genome; should not include any file extensions
5) email: Optional; input your email and Westgrid will supposedly notify you when the job has finished

Outputs:
1) PBS files; PBS for each file to be analyzed (FtbA) will have the name <FtbA>.pbs
2) Script will make directories to save the analysis results in if they are not already made

'''

from os.path import join, isfile, isdir, basename, dirname
from os import mkdir, listdir
import sys, argparse


def Directory_Maker(directory_list):
	for folder in directory_list:
		if not isdir(folder):
			print("Making the directory %s" %folder)
			mkdir(folder)


def BT_ST_PBS_Maker(results_dir, data_file, reference_genome, email):
	#Define and make directories in which to save results
	trimmed_dir = join(results_dir,'trimmed_seqs')
	bowtie_dir = join(results_dir, 'bowtie_alignments')
	all_reads_dir = join(bowtie_dir, 'all_mapped_reads')
	unique_reads_dir = join(bowtie_dir, 'unique_map_reads')
	multi_map_reads_dir = join(bowtie_dir, 'multi_map_reads')
	bam_list_dir = join(results_dir, 'bam_file_lists')
	mpileup_dir = join(results_dir, 'mpileup_results')
	directory_list = [results_dir, trimmed_dir, bowtie_dir, all_reads_dir, unique_reads_dir,
	 multi_map_reads_dir, bam_list_dir, mpileup_dir]
	Directory_Maker(directory_list)

	#Define file names for each result/step
	sample_name = basename(data_file).replace('.fq.gz', '')
	trim_seq = join(dirname(reference_genome), 'TruSeq3-SE.fa')
	trimmed_file = join(trimmed_dir, sample_name + '_trimmed.gz')
	bowtie_file_out = join(bowtie_dir, sample_name + '.sam')
	init_bam_file = join(all_reads_dir, sample_name + '.bam')
	sorted_bam_file = join(all_reads_dir, sample_name + '_sorted.bam')
	rmdup_bam_file = join(all_reads_dir, sample_name + '_rmdup.bam')
	unique_bam_file = join(unique_reads_dir, sample_name + '_unique_read.bam')
	multi_map_bam_file = join(multi_map_reads_dir, sample_name + '_multi_map_read.bam')
	full_bam_list = join(bam_list_dir, "full_bam_list.txt")
	unique_bam_list = join(bam_list_dir, "unique_bam_list.txt")
	multi_bam_list = join(bam_list_dir, "multi_bam_list.txt")

	#Assemble PBS file piecewise for easier variable referencing
	pbs_part1 = """####
##
## Resource Allocation Section
##
####

## Set the shell
#! /bin/bash
#PBS -S /bin/bash

##Main request settings:
#PBS	-l	nodes=1:ppn=2
#PBS	-l	walltime=24:00:00
#PBS	-l	mem=4gb

##Email settings
#PBS -M %s
#PBS -m abe

####
##
## The Script Portion
##
####

#Change to appropriate directory
cd $PBS_O_WORKDIR

## Initialize required programs
module load java
module load application/bowtie2/2.2.3
module load application/samtools/1.3.1
JAR=/global/software/trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar
""" %email
	pbs_part2 = "\n##Run trimmomatic\n" + \
	"java -jar $JAR SE -threads $PBS_NUM_PPN -phred33 %s %s " %(data_file, trimmed_file) + \
	"ILLUMINACLIP:%s:2:30:6 HEADCROP:20 SLIDINGWINDOW:4:15 MINLEN:36 CROP:170 \n" %(trim_seq)  
	pbs_part3 = "\n## Run Bowtie2\n" + \
	"bowtie2 --end-to-end --sensitive -p $PBS_NUM_PPN --non-deterministic " + \
	"-x %s -U %s -S %s \n" %(reference_genome, trimmed_file, bowtie_file_out)
	pbs_part4 = "\n## Clean-up and format data with Samtools\n" + \
	"samtools view -bS %s > %s \n" %(bowtie_file_out, init_bam_file) + \
	"samtools sort %s -o %s \n" %(init_bam_file, sorted_bam_file) + \
	"samtools rmdup -s %s %s \n" %(sorted_bam_file, rmdup_bam_file) + \
	"rm %s \n" %init_bam_file + \
	"mv %s %s \n" %(rmdup_bam_file, init_bam_file)
	pbs_part5 = "\n## Break results into unique and multi-mapped reads\n" + \
	"samtools view -b -q 10 %s > %s \n" %(init_bam_file, unique_bam_file) + \
	"samtools view -b -q 10 -U %s > %s \n" %(init_bam_file, multi_map_bam_file)
	pbs_part6 = "\n## Index all BAM files for referencing/viewing\n" + \
	"samtools index %s \n" %init_bam_file + \
	"samtools index %s \n" %unique_bam_file + \
	"samtools index %s \n" %multi_map_bam_file
	pbs_part7 = "\n## Add files to bam file list (for use with mpileup)\n" + \
	"echo \"%s\" >> %s \n" %(init_bam_file, full_bam_list) + \
	"echo \"%s\" >> %s \n" %(unique_bam_file, unique_bam_list) + \
	"echo \"%s\" >> %s \n" %(multi_map_bam_file, multi_bam_list)
	pbs_part8 = "\n## Clean up the created files\n" + \
	"rm -f %s \n" %bowtie_file_out + \
	"rm -f %s \n" %sorted_bam_file
	
	total = pbs_part1 + pbs_part2 + pbs_part3 + pbs_part4 + pbs_part5 + pbs_part6 + \
	pbs_part7 + pbs_part8
	return(total)

def File_Writer(text, file_name):
	with open(file_name, 'w') as f:
		f.write(text)

def __Main__(data_dir, output_dir, results_dir, referene_genome, email = 'blank@email.com'):
	file_list = [join(data_dir, f) for f in listdir(data_dir) if isfile(join(data_dir,f))]
	Directory_Maker([output_dir])
	for file in file_list:
		out_file = join(output_dir ,basename(file).replace(".fa.gz", "") + ".pbs")
		text_string = BT_ST_PBS_Maker(results_dir, file, referene_genome, email)
		File_Writer(text_string, out_file)

def Arg_Parse_Init():
	parser = argparse.ArgumentParser(description='Write PBS scripts for the bowtie and samtools alignment process')
	parser.add_argument('--data', '-d', help="path to directory containing data", required=True)
	parser.add_argument('--results_dir', '-r', help="path to directory where BT_SP pipeline results will be saved", required=True)
	parser.add_argument('--pbs_dir', '-o', help="path where PBS script files will be saved; default is current directory", default='./')
	parser.add_argument('--ref_genome', '-g', help="path to the reference genome, do not include any file extensions", required = True)
	parser.add_argument('--email', '-e', help="send process updates to this email address")
	return parser

def Action_Tree(args):
	if not args.email:
		__Main__(args.data, args.pbs_dir, args.results_dir, args.ref_genome)
	else:
		__Main__(args.data, args.pbs_dir, args.results_dir, args.ref_genome, args.email)

args = Arg_Parse_Init().parse_args()
Action_Tree(args)

