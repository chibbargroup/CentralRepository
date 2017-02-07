'''
Transpose.py
Written by J. Hayes 2.6.2017

Purpose:
Reads data taken from the VESPERS beamline and transposes it for processing in PyMCA.

Input:
Program is self-contained, but will ask for the directory containing the data to transpose
and an output directory to save the tranposed data; note that the output directory cannot be the input directory

Output:
Saves the tranposed data as .csv files in the specified output directory

'''

from os.path import join, isfile, isdir
from os import listdir, mkdir
import pandas as pd
import re
def Boiler_Plate():
	print('*'*61)
	print('*'+' '*59 +'*')
	print('*' + ' '*25 + 'TRANSPOSE' + ' '*25 + '*')
	print('*'+ ' '*int((59-len('v 0.1, (C) Feb 2017'))/2) + 'v 0.1, (C) Feb 2017' + ' '*int((59-len('v 0.1, (C) Feb 2017'))/2) + '*')
	print('*'+' '*59 +'*')
	print('*'*61)

def Get_Data_Directory():
	good_input = False
	while not good_input:
		directory = input("Please input the directory containing the files to transpose: ")
		if isdir(directory):
			good_input = True
		else:
			print("Hmm...that directory doesn't seemt to exist, please try again")
	return directory

def Get_Output_Directory(data_dir):
	good_input = False
	while not good_input:
		directory = input("Please input the directory where you'd like me to save the transposed files: ")
		if directory != data_dir:
			if not isdir(directory):
				print("Creating directory %s" %directory)
				mkdir(directory)
			good_input = True
		else:
			print("Sorry, the output directory has to be different from the input directory")
	return directory

def Transpose_Data(data_dir, output_dir):
	file_list = [f for f in listdir(data_dir) if isfile(join(data_dir, f))]
	for file in file_list:
		df = pd.read_csv(join(data_dir, file), delimiter = '\t', index_col = False, header = None)
		df = df.transpose()
		new_file_name = join(output_dir, re.sub('_.+_', '_', file))
		print("Writing...%s" %new_file_name)
		df.to_csv(new_file_name, index = False, header = False)

def __Main__():
	data_dir = Get_Data_Directory()
	output_dir = Get_Output_Directory(data_dir)
	Transpose_Data(data_dir, output_dir)

Boiler_Plate()
__Main__()