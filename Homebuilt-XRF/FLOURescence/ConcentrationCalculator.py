'''
Concentration Calculator
Written by J. Hayes, Last Editted 1/27/2017

Purpose: This script uses the calibration curves to calculate the concentration of each analyte
for each sample. The script also calculates the error associated with each measurement. 

Input:
1) raw_sample_data_file - contains the combined, unaveraged sample data
2) avg_sample_data_file - contains the combined and averaged sample data
3) calibration_curve_file - the csv file containing the calibration curve info
4) output_dir - the directory where the results will be saved

Output:
The script generates two csv files which are saved in the specified output_dir
1) sample_concentrations.csv: contains the un-averaged sample concentrations
2) avg_results.csv: contains the averaged sample concentrations and the corresponding error in each
measurement

'''

import pandas as pd 
import numpy as np 
from os import listdir
from os.path import isfile, join
from collections import Counter

#Uses the calibration curve data to calculate the concentration of each sample
def Concentration_Calculator(data, cal_curves):
	for element in data.index:
		data.ix[element] = (data.ix[element]-cal_curves[element]['intercept'])/cal_curves[element]['slope']
	return(data)

#Counts the number of times each sample was replicated
def Replicate_Counter(data):
	sample_list = []
	for column in data:
		column = column.split('.')[0]
		sample_list.append(column)
	replicate_count = Counter(sample_list)
	return replicate_count

#Calculates the error for each measurement and then adds it to the main dataframe
def Error_Analysis(data, cal_curves, replicate_count):
	sample_error = {}
	for column in data:
		single_error = {}
		for element in data[column].index:
			x_avg = data.ix[element,column]
			std_dev_y = cal_curves.ix['std_dev_y', element]
			xi = cal_curves.ix['xi', element]
			xi_2 = cal_curves.ix['xi_2', element]
			m = abs(cal_curves.ix['slope', element])
			replicate = replicate_count[column]
			n = cal_curves.ix['n', element]
			d = n*xi_2 - xi**2
			error = std_dev_y/m*(1/replicate + (x_avg*n)/d + xi_2/d + 2*x_avg*xi/d)**0.5
			single_error[element] = error
		sample_error[column + '-Error'] = single_error
	error_results = pd.DataFrame.from_dict(sample_error)
	final_results = pd.concat([data, error_results], axis = 1)
	final_results = final_results.sort_index(axis=1)

#Runs the script
def Concentration_Analysis(raw_sample_data_file, avg_sample_data_file, calibration_curve_file, output_dir):
	#Open relevant files as dataframes
	cal_curves = pd.read_csv(calibration_curve_file, index_col = 0)
	replicate_count = Replicate_Counter(raw_sample_data_file)
	raw_data = pd.read_csv(raw_sample_data_file, index_col = 0)
	avg_data = pd.read_csv(avg_sample_data_file, index_col = 0)
	
	#Analyze raw data and save as CSV file (asked for by biologists)
	raw_concs = Concentration_Calculator(raw_data, cal_curves)
	replicate_count = Replicate_Counter(raw_data)
	raw_concs.to_csv(join(output_dir, "sample_concentrations.csv"))
	
	#Analyze the avg data points
	avg_data = Concentration_Calculator(avg_data, cal_curves)
	results = Error_Analysis(avg_data, cal_curves, replicate_count)
	results.to_csv(join(output_dir, "avg_results.csv"))
	return results
