'''
Concentration and Error Analysis Script
Written by J.Hayes; Last Edited Oct. 26, 2016

Purpose: Calculates the analyte concentration and performs error analysis based on the results of CalibrationCurves.py
Usage: Requires data in CSV file, npy array with concentrations, linear regression stats, and combined calibration data

Input: datafiledir = the directory containing all the data that was output by PlateData.py
	   calinfofile = the csv file containing the linear regression parameters output by CalibrationCurves.py
	   caldata = the csv file containing the combined calibration peak areas
	   calvalues = the npy array that lists the concentration of analyte for each calibration sample
	   outputdirectory = the place where the output should be saved
Output: Regression_Stats.txt = text file containing the tabulated regression stats (i.e., std_y, sumx, D, etc)
        Analysis_Results.csv = csv file that tabulates the calculated element concentration and the error associated with that value


'''

import numpy as np
import pandas as pd
from os import listdir, walk, mkdir
from os.path import isfile, join, exists, abspath, isdir
from collections import Counter

#Opens the data files and combines the data. Then removes the superfluous data sets
def DataCombiner(datafiledir):
	combdata = pd.DataFrame()
	
	#opens all the data-containing files and puts them into one dataframe
	all_files = all_files = [f for f in listdir(datafiledir) if isfile(join(datafiledir, f))]
	for f in all_files:
		if "peakfitnorm" in f:
			df = pd.read_csv(join(datafiledir,f), index_col = 0)
			combdata = pd.concat([combdata,df], axis = 1)
	
	#removes the spectra collected on empty spaces and calibration samples
	cols = set(combdata.columns)
	for i in cols:
		if 'Cal ' in i or 'x' in i or 'empty' in i or 'Empty' in i or 'X' in i:
			del combdata[i]

	#try to rename columns		
	holder = []
	for col in combdata.columns:
		col = str(col)
		#if the column isn't in the holder (i.e., is currently unique), add it to the lsit
		if col not in holder:
			holder.append(col)
		#if holder isn't unique, then make it unique
		else:
			found = False
			for i in range(0,9):
				if len(col) < 3:
					col = col + '.1'
				elif col[-2] != '.':
					col = col + '.1'
				elif col[:-1] + str(i + 1) not in holder:
					holder.append(col[:-1] + str(i +1))
					found = True
				
				if found:
					break

	combdata.columns = (holder)

	combdata = combdata.transpose()

	del combdata['Io']
	
	combdata.to_csv('test.csv')
	return(combdata)

			
def ConcCalc(data,calinfofile):
	combdata = pd.DataFrame()
	calstats = pd.read_csv(calinfofile, index_col = 0)
	for column in data:
		conccol = (data[column]-calstats.ix['intercept',column])/calstats.ix['slope',column]
		combdata = combdata.append(conccol)
		combdata.to_csv('test.csv')
	return(combdata)

#Determines which samples are repeats, and averages the values
#Note that method assumes no more than 10 repetitions per sample
def Data_Average(data):
	sample_ids = list(data.columns.values)
	renamed = []
	for i in sample_ids:
		if len(i) < 3:
			i = '__' + i
		
		renamed.append(i)

	bases = [i for i in renamed if i[-2] != '.']
	repeats = [i for i in renamed if i[:-2] in sample_ids]
	combdata = pd.DataFrame()

	#Note: method assumes no more than 10 reptitions per sample
	for i in bases:
		j = 1
		for y in repeats:
			if y[:-2] == i:
				j += 1
				total = data.ix[:,i]+data.ix[:,y]
		avg = total/j
		avg.index.rename(i)
		
		holder = pd.DataFrame(avg, columns = [i])
		counter = pd.DataFrame(data = [j], columns = [i], index = ['N'])
		holder = pd.concat([holder,counter], axis = 0)
		combdata = pd.concat([combdata, holder], axis = 1)
	return(combdata)

#Calculates the regression statistics needed for the error analysis; also saves the error analysis stats in a text file
#Items returned 1) dic, a dictionary containing N, D, sumx, and sumxsquare; and 
#2) df2, a data frame containing the LOD, LOQ, and std_y values
def RegressStats(calinfofile,caldata,calvalues,dir_out):
	#Read data in
	Reg = pd.read_csv(calinfofile, index_col = 0)
	Data = pd.read_csv(caldata, index_col = 0)

	#Calculate N
	N = len(Data.columns.values)

	#Calculates LOD and LOQ for each element
	LOD = Reg.ix['std_err']*3.3/Reg.ix['intercept']
	LOQ = Reg.ix['std_err']*10/Reg.ix['intercept']
	
	#Calculates sum(x) and sum((x)^2), D
	sumxsquare = 0
	sumx = 0
	column_list = [i[:5] for i in (Data.columns.values)]
	cal_count = Counter(column_list)
	cal_concs = np.load(calvalues)
	i = 1
	for conc in cal_concs:
		sumx += conc * cal_count['Cal ' + str(i)]
		sumxsquare += conc**2 * cal_count['Cal ' + str(i)]	
		i += 1
	D = N*sumxsquare - sumx**2

	#Calculate di2 and std_y
	Reg_Transpose = Reg.transpose()
	df = pd.DataFrame()
	for col in Data:
		conc = cal_concs[int(col[4])-1]
		holder = (Data[col]-Reg_Transpose['slope']*conc-Reg_Transpose['intercept'])**2
		df = pd.concat([df,holder], axis = 1) ##This is the step where the order of the indexes gets messed up; fix here?
	df['sumdi2'] = df.sum(axis = 1)
	sumdi2 = df['sumdi2']
	std_y = (sumdi2/(N-2))**0.5

	df2 = pd.concat([std_y, LOD, LOQ], axis = 1)
	df2.columns = ['std_y', 'LOD', 'LOQ']

	#write regression stats to file
	with open(join(dir_out,'Regression_Stats.txt'), 'w') as f:
		f.write('N = ' + str(N) + '\n')
		f.write('D = ' + str(D) + '\n')
		f.write('sumx = ' + str(sumx) + '\n')
		f.write('sumxsquare = ' + str(sumx) + '\n')
		f.write('std_y = \n' + str(std_y.rename()) + '\n')
		f.write('LOD = \n' + str(LOD) + '\n')
		f.write('LOQ = \n' + str(LOQ) + '\n')

	dic = {'N': N, 'D': D, 'sumx': sumx, 'sumxsquare': sumxsquare}
	return dic, df2

def Error_Analysis(data, reg_int, reg_arrays, reg_stats, file_out): 
	#open the reg_stats and transpose into easily usable data array; this step could likely be avoided if we add an open data function
	reg_stats = pd.read_csv(reg_stats, index_col = 0).transpose()
	
	#This block of code largely made to deal with the mis-ordering of reg_arrays; if that problem could be fixed
	#this code is likely not needed
	reg_stats = reg_stats.sort_index()
	N = data.ix['N', :]
	data_no_n = data.drop(data.index[-1]).sort_index()
	#Set holder, data labels
	error_df = pd.DataFrame()
	data_labels = data.columns.values

	#Calculates the error for each value
	for j in range(0, len(data.columns.values)):
		error_out_of_root = reg_arrays['std_y']/reg_stats['slope']
		error_under_root = 1/data.ix['N', j] + data_no_n.ix[:,j]**2*reg_int['N']/reg_int['D']+reg_int['sumxsquare']/reg_int['D']+2*data_no_n.ix[:,j]*reg_int['sumx']/reg_int['D']
		error = error_out_of_root*error_under_root**0.5
		error = error.rename(data_labels[j])
		error_df = pd.concat([error_df, error], axis = 1)

	df = pd.concat([data_no_n, error_df], axis = 1)
	
	#Rename duplicate columns and sort so data and error columns are paired
	cols=pd.Series(df.columns)
	for dup in df.columns.get_duplicates():
			cols[df.columns.get_loc(dup)]=[dup+'.'+str(d_idx) if d_idx!=0 else dup for d_idx in range(df.columns.get_loc(dup).sum())]
	df.columns=cols
	df = df.reindex_axis(sorted(df.columns), axis=1)

	#Determine what values are below LOD, LOQ; this part is currently busted
	for j in range(0, len(df.columns.values), 2):
		i = 0
		for x in df.ix[:, j]:
			if x < reg_arrays.ix[i,'LOD']:
				df.ix[i,j] = 'Below LOD'
				df.ix[i,j+1] = 'Below LOD'
			elif x < reg_arrays.ix[i, 'LOQ']:
				df.ix[i,j] = 'Below LOQ'
				df.ix[i,j+1] = 'Below LOQ'
			i += 1

	#Rename the error columns 
	new_names = []
	for col in df.columns.values:
		if '.1' not in col:
			new_names.append(col)
		else:
			new_names.append('Error')
	df.columns = new_names
	
	df.to_csv(join(file_out,'Analysis_Results.csv'))

def Analyzer(data_dir, cal_file, cal_comb, cal_concs, file_out):
	print('Starting Concentration Calculator')
	print('Renticulating Splines...')
	combdata = DataCombiner(data_dir)
	conc_data = ConcCalc(combdata, cal_file)
	avg_data = Data_Average(conc_data)
	reg_stats, loq_lod = RegressStats(cal_file, cal_comb, cal_concs, file_out)
	Error_Analysis(avg_data, reg_stats, loq_lod, cal_file, file_out)

#For testing purposes

# datafiledir = 'C:/Users/John/Desktop/tester/Norm'
# calinfofile = 'C:/Users/John/Desktop/tester/Calibration_Curves/CalibrationCurvesOut.csv'
# calcomb = 'C:/Users/John/Desktop/tester/Calibration_Curves/CombinedCals.csv'
# cal_concs = 'C:/Users/John/Desktop/tester/Calibration_Curves/concentrations.npy'
# file_out = 'C:/Users/John/Desktop/tester/'
# Analyzer(datafiledir, calinfofile, calcomb, cal_concs, file_out)