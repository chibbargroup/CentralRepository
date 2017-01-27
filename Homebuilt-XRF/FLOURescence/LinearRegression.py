"""
Linear_Regession.py
Written by D. Sanche and J. Hayes, Last Edit: 1/27/17

Purpose: Calculate the calibration curves via linear regression of the calibration data

Input:
1) cal_data_file: data file containing the calibration data
2) output_dir: where to save the output results

Note that the script will look for the cal_concentrations file in the output directory; if it finds
this file it will ask if the user would like to use these values for the calibration curve analysis

Output:
Saves the regression data to file calibration_curves.csv file located in the output directory
This file contains the following: slope, intercept, r-value, p-value, std error, # of pts, 
LOD, LOQ, xi^2, and xi

Also saves the calibration values to the file cal_concentrations in the specified output directory

"""
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
from scipy import stats

#Gets the calibration concentrations as input from user; saves inputs as csv
def Get_Cal_Concs_From_User(data_file, output_dir):
	df = pd.read_csv(data_file, index_col = 0)
	concentrations = {}
	for cal in df.columns: 	
		cal_conc = Get_User_Input(cal)
		concentrations[cal] = cal_conc
	
	concentrations = pd.Series(concentrations).rename('Concentrations')
	concentrations.to_csv(join(output_dir, 'cal_concentrations.csv'))
	df = df.append(concentrations)
	return(df)

#Gets concentration inputs from user; ensures user can only input numbers
def Get_User_Input(calibration_number):
	is_number = False
	while not is_number:
		number = input("Please input the concentration of %s: " %calibration_number)
		try:
			float(number)
			is_number = True
		except ValueError:
			is_number = False
			print("I'm sorry Dave, I can't do that. It seems you didn't input a number...please concentrate")
	return(float(number))

#Gets the calibration concentrations from the previously saved csv file
def Get_Cal_Concs_From_File(data_file, conc_file):
	df = pd.read_csv(data_file, index_col = 0)
	concentrations = pd.read_csv(conc_file, names = ['', 'Concentrations'], index_col = 0).transpose()
	df = df.append(concentrations.ix['Concentrations'])
	return(df)

#Adds calibration concentrations to dataframe containing calibration data
def Add_Cal_Concs(data_file, output_dir):
	#Initialize some parameters/response lists
	yea_response = ['Y', 'y', 'yes', 'Yes']
	nay_response = ['N', 'n', 'No', 'no']
	concentration_file = join(output_dir, 'cal_concentrations.csv')

	#Run the logic flow
	readable_input = False
	if isfile(concentration_file):
		while not readable_input:
			use_previous = input("I see that there is already a concentration file, would you like me to use it? [Y/N]: ")
			if use_previous in yea_response:
				df = Get_Cal_Concs_From_File(data_file, concentration_file)
				readable_input = True
			elif use_previous in nay_response:
				df = Get_Cal_Concs_From_User(data_file, output_dir)
				readable_input = True
			else:
				print("I'm sorry Dave, I can't do that. Please try answering again.")
	else:
		df = Get_Cal_Concs_From_User(data_file, output_dir)
	return(df)

#The hear to the program; runs the linear regression analysis for each element and returns
#the results as a dataframe
def Linear_Regression_Calculator(cal_data):
	cal_curve_list = []
	concentrations = np.array(cal_data.ix['Concentrations'])
	for element in cal_data.index:
		if element != "Concentrations":
			data = np.array(cal_data.ix[element])
			slope, intercept, r_value, p_value, std_error = stats.linregress(concentrations, data)
			n, lod, loq, xi_2, xi = Regression_Stats(concentrations, std_error, intercept)
			std_dev_y = Std_Dev_Y_Calculator(data, concentrations, slope, intercept)
			reg_values = [slope, intercept, r_value, p_value, std_error, n, lod, loq, std_dev_y, xi, xi_2]
			reg_labels = ["slope", "intercept", "r_value", "p_value", "std_err", "n", "LOD", "LOQ", "std_dev_y", "xi", "xi_2"]
			cal_curve = pd.DataFrame(reg_values, columns=[element], index=reg_labels)
			cal_curve_list += [cal_curve]
	cal_curve_results = pd.concat(cal_curve_list, axis = 1)
	return cal_curve_results

#Separate script to calculate various extra statistics about the regression;
#runs as part of the Linear_Regression_Calculator module
def Regression_Stats(concentrations, std_error, intercept):
	n = len(concentrations)
	lod = 3.3*std_error/intercept
	loq = 10*std_error/intercept
	xi_2 = 0
	xi = 0
	for conc in concentrations:
		xi += conc
		xi_2 += conc**2
	return n, lod, loq, xi_2, xi

#Calculate the standard deviation in the y values
def Std_Dev_Y_Calculator(cal_data, concentrations, slope, intercept):
	di_2 = 0
	for i in range(0, len(cal_data)):
		di = cal_data[i] - slope*concentrations[i] + intercept
		di_2 += di**2
	std_dev_y = (di_2/(len(concentrations)-2))**0.5
	return std_dev_y

#Run the script
def Linear_Regression_Script(cal_data_file, output_dir):
	data = Add_Cal_Concs(data_file, output_dir)
	cal_curves = Linear_Regression_Calculator(data)
	cal_curves.to_csv(join(output_dir, "calibration_curves.csv"))

Linear_Regression_Script(data_file, output_dir)
