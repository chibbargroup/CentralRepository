'''
NormalizeData 
Written by J.Hayes
Last editted 10/24/2016

Program adds header labels and normalizes fit data processed from FormatPlate.py
Requires 3 inputs:
1) Directory with peak fit data
2) Directory with Io files
3) File containing header information for all plates (in CSV format); see program documentation for more into on the format of this file

'''
import numpy as np
import pandas as pd
from os import listdir, walk, mkdir
from os.path import isfile, join, exists, abspath, isdir
import sys

#Takes a single data file and the header file and adds sample labels based on the file name; returns the data in a dataframe
def Replace_Header(datafile, headerfile):
	data = pd.read_csv(datafile)
	header = pd.read_csv(headerfile)
	lb_array = header['Plate']
	header = header.set_index(lb_array)
	for n in header.ix[:,0]:	
		if n in datafile:
			holder = list(header.ix[n])
			data.columns = holder
	return(data)

	balls = list(newfile.ix[0,:])
	newfile.columns = balls

#Gets Io data from Io file
def Io_Strip(IoFile):
	#Open file containing Io data and extract the Io values to a list
	print(IoFile)
	file = open(IoFile, 'r')
	IoData = []
	for line in file:
		if '#' not in line:
			linelist = [float(i) for i in line.strip().split()] 
			if len(linelist) > 5:
				IoData.append(linelist[4])
	return(IoData)

#Reads the data and the Io and puts both into sets into a single dataframe
def Add_Io_to_Peakfiles(IoFileDir, Peak_data, Peak_file_Name):
	#Generate list of Io Files
	all_files = [f for f in listdir(IoFileDir) if isfile(join(IoFileDir, f))]
	for f in all_files:
		print(Peak_file_Name.strip('_Ge13El_1.dat_FITDIR_peakfit.csv'))
		if Peak_file_Name.strip('_Ge13El_1.dat_FITDIR_peakfit.csv') in f: 
			Io_values = Io_Strip(join(IoFileDir,f))
			Io_values.insert(0, 'Io')
			df = pd.DataFrame(Io_values)
			df = df.transpose()
			df.columns = Peak_data.columns.values
			Peak_data = Peak_data.append(df, ignore_index = True)
			return(Peak_data)

#Element selector reads a dataframe and selects only the elements to analyze for export; returns data in a dataframe
def Element_Selector(CombData, is_calcium):
	if is_calcium:
		Element_List = ['Ca', 'Io']
	else:
		Element_List = ['Mn','Cu','Co','Ni','Fe','Zn','Mo','Se','Rb','Cd','Io']	
	i = 0
	for el in CombData.ix[:,0]:
		if el not in Element_List:
			CombData = CombData.drop(i)
		i += 1
	CombData = CombData.reset_index()
	del CombData['index']
	return CombData

#Normalize_Area reads a data frame and then normalizes the data; returning the data in a dataframe
def Normalize_Areas(RawPeakAreas):
	IoRow = len(RawPeakAreas.index)-1
	El_Labels = np.asarray(RawPeakAreas.ix[:,0])
	del RawPeakAreas[RawPeakAreas.columns[0]]

	Normalized = RawPeakAreas/RawPeakAreas.values[IoRow,:]

	Normalized.insert(0, 'Element Label', El_Labels)
	return(Normalized)

#Opens all the peakfit files in a directory, adds headers to the files, and then saves the result
def Replace_Header_Batch(datafiledir, iofiledir, headerfile, output_dir, calcium):
	savedir = join(output_dir, 'Norm')

	if not exists(savedir):
		mkdir(savedir)
	
	all_files = [f for f in listdir(datafiledir) if isfile(join(datafiledir, f))]
	for name in all_files:
		if "peakfit" in name and "peakfitnorm" not in name:
			f = join(datafiledir,name)			
			newfile = Replace_Header(f, headerfile)
			newfile = Add_Io_to_Peakfiles(iofiledir, newfile, name)
			newfile = Element_Selector(newfile, calcium)
			newfile = Normalize_Areas(newfile)
			newfile.to_csv(join(savedir, name.strip('.csv')) + 'norm.csv', index = False)
			print('Writing...' + join(savedir, name.strip('.csv')) + 'norm.csv')

