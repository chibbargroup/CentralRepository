'''

'''

import numpy as np
import pandas as pd
from os import listdir, walk, mkdir
from os.path import isfile, join, exists, abspath, isdir
from collections import Counter

#Opens and combines the data. Then removes the superfluous data sets
def Data_Combiner(norm_file_dir):
	data = pd.DataFrame()
	
	#opens all the data-containing files and puts them into one dataframe
	all_files = all_files = [f for f in listdir(norm_file_dir) if isfile(join(norm_file_dir, f))]
	for f in all_files:
		df = pd.read_csv(join(datafiledir,f), index_col = 0)
		data = pd.concat([data,df], axis = 1)
	return data 

#removes the spectra collected on empty spaces and calibration samples
def Column_Cleaner(data):
	empty_aliases = ['empty', 'Empty', 'x', 'X', 'blank', 'Blank'] #Can be modified as necessary
	for column in set(data.columns): #Note: Use set because del data[column] will remove all instances of that column value
		column_read = column.split('.')[0]	
		if column_read in empty_aliases:
			del data[column]
		elif "Cal " in column:
			del data[column]
	return data
	
def Replicate_Renamer(data):
	column_list = list(data.columns)
	column_counts = Counter(column_list)
	holder = []
	test = list(set(column_list))
	for column in column_list:
		rename = False
		i = 1
		if column in test:
			rename = True
			while rename:
				if '.' in column:
					column = column.split('.')[0]
				new_name = column + '.%s' %i
				if new_name not in test:
					holder.append(new_name)
					test.append(new_name)
					rename = False
				else:
					i += 1
	print(holder)

	'''
	for column in column_list:
		if column_counts[column] > 1:
			print(column)
			print(column_counts[column])
			rename = True
			i = 1
			while rename:
				if '.' in column:
					print(type(column))
					j = column.split('.')
					
				new_name = column + ".%s" %i
				if new_name in column_counts.keys():
					i += 1
				else:
					rename = False
			#column_list[j] = new_name
			column_counts = Counter(column_list)

	keys = list(test.keys())
	for key in test:
		if test[key] > 1:
			holder.append(key)
	
	for thing in holder:
		i = 1
		while thing in keys:
			thing = thing + ".%s" %i
			i += 1
		keys.append(thing)
	'''


datafiledir = 'C:/Users/John/Desktop/Test/FIT'
calinfofile = 'C:/Users/John/Desktop/Test/Calibration_Curves/CalibrationCurvesOut.csv'
calcomb = 'C:/Users/John/Desktop/Test/Calibration_Curves/CombinedCals.csv'
cal_concs = 'C:/Users/John/Desktop/Test/Calibration_Curves/concentrations.npy'
file_out = 'C:/Users/John/Desktop/Test/'

df = Data_Combiner(datafiledir)
df = Column_Cleaner(df)
Replicate_Renamer(df)