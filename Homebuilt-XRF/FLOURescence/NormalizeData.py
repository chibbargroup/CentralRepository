'''
NormalizeData_v2
Written by J.Hayes
Last editted 1/24/2017

Program reads through peakfit files, adds spectra labels and io data, and then returns the
normalized peak areas for the elements of interest

Program Inputs:
1) peakfit file directory
2) header file
3) io file directory
4) output directory 
5) calcium analysis boolean (true or false)

Program Outputs:
None (saves normalized data in Norm folder created in output directory)

'''

from os import listdir, mkdir
from os.path import isfile, join, isdir, basename
import numpy as np
import pandas as pd

#Replaces the headers so that each set of peak fits is labelled with the sample name
def Replace_Header(data_file, header_file):
	data = pd.read_csv(data_file, index_col = 0)
	headers = pd.read_csv(header_file, index_col = 0)
	for plate in headers.index:
		if plate in data_file:
			data.columns = list(headers.ix[plate])
	return data

#Parse through all elements investigated and return only elements of interest (dataframe)
def Element_Parse(data, ca_analysis):
	if ca_analysis:
		element_list = ['Ca']
	else:
		#Can change list as needed
		element_list = ['Mn','Cu','Co','Ni','Fe','Zn','Mo','Se','Rb','Cd'] 

	for element in data.index:
		if element not in element_list:
			data = data.drop(element)
	return data

#Searches through list of Io files and picks appropriate Io file
def Io_Grabber(data_file, io_dir):
	io_files = [f for f in listdir(io_dir) if isfile(join(io_dir, f))]
	data_file = basename(data_file).replace("_Ge13El_1_peakfit.csv", "") + "_"
	
	for file in io_files:
		if data_file in file:
			io_data = Io_File_Parse(join(io_dir, file))
			break
	return io_data

#Parses through the Io file and picks out the Io values, returns them as a data frame
def Io_File_Parse(io_file):
	io_data = []
	with open(io_file, 'r') as f:
		for line in f:
			if "#" not in line and line != '\n':
				line = line.split('\t')
				io_value = float(line[4]) #Number here is Io column; may change
				io_data.append(io_value)
	df = pd.DataFrame({'Io': io_data})
	df = df.transpose()
	return df

#Normalizes the data by dividing each peak area by the corresponding Io value
def Normalize_Areas(data, io):
	io.columns = data.columns.values
	data = data.append(io)
	data = data/data.ix['Io']
	return data

def Data_Process(data_file, header_file, io_dir, ca_analysis):
	df = Replace_Header(data_file, header_file)
	df = Element_Parse(df, ca_analysis)
	io = Io_Grabber(data_file, io_dir)
	df = Normalize_Areas(df, io)
	return df

def Batch_Process(peakfit_dir, header_file, io_dir, output_dir, ca_analysis):
	#Make sure save directory exists, if not, make it
	save_dir = join(output_dir, 'Norm')
	if not isdir(save_dir):
		mkdir(save_dir)

	file_list = [f for f in listdir(peakfit_dir) if isfile(join(peakfit_dir, f))]
	
	for file in file_list:
		save_name = file.replace('.csv', '') + '_norm.csv'
		save_name = join(save_dir, save_name)
		file = join(peakfit_dir, file)
		peak_fit = Data_Process(file, header_file, io_dir, ca_analysis)
		peak_fit.to_csv(save_name)

