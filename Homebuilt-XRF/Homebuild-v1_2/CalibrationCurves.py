"""
Created by Daniel Sanche on October 24, 2016
Purpose: finds the "Cal" columns of every .csv file in a directory, and performs linear regressions for each element in the files
Usage: python CalibrationCurves.py input_files_dir output_files_dir concentration_file.npy
Inputs:
    input_files_dir:        the path to the directory that contains the input .csv files
                            if left blank, the script's parent directory is used
    output_files_dir:       the path to the folder where results will be saved
                            if left blank, the script's parent directory is used
    concentration_file.npy  represents the concentration for each Cal value
                            if left blank, the program will ask the user to input all values
                            if the script has been run before, you can use the concentration file generated last time
Outputs:
    CalibrationCurvesOut.csv     contains the results of the linear regression for each element
    CombinedCals.csv                a single .csv file that represents the "Cal N" columns of all .csv files in the input_files_dir
    concentrations.npy              a saved version of the the user input, which can be input to the next run to save time
"""

import sys
import os
import  pandas as pd
import numpy as np
from scipy import stats

#returns whether a string can be interpreted as a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#extracts the Cal columns from a .csv file
def getCalsFromFile(file_path, dropRows=["Io"]):
    num_cals = 0
    df = pd.read_csv(file_path, index_col=0)
    for column in df.columns.values:
        if "Cal" not in column:
            df = df.drop(column, 1)
        else:
            numStr = column.replace("Cal ", "")
            num_cals = max(num_cals, int(numStr))
    #remove unwanted rows
    df = df.drop(dropRows)
    return df, num_cals

#extracts a list of dataframes with only the Cal columns from a directory of .csv files
#also returns the max number of Cal columns in any .csv file
def getCalsFromDir(dir="./"):
    num_cals = 0
    dataframeList = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if ".csv" in file:
                this_df, this_num = getCalsFromFile(os.path.join(root,file))
                num_cals = max(num_cals, this_num)
                dataframeList += [this_df]
    return dataframeList, num_cals

#asks the user to input concentrations for each Cal value
def inputConcentrations(num_cals):
    results = []
    for i in range(num_cals):
        numberInput = False
        while not numberInput:
            concentaration = input("Enter Concentration for Cal " + str(i+1) + ":   ")
            numberInput = is_number(concentaration)
            if not numberInput:
                print ("Error: " + str(concentaration) + " is not a number")
        results += [float(concentaration)]
    return results

#performs a linear regression from a list of dataframes containing Cal columns
#returns a dataframe of the results for each element
def performLinearRegression(dataframeList, concentrations, rWarnVal = 0.9):
	combined = pd.concat(dataframeList, axis=1)
	resultList = []
	for elem in combined.index:
		xVals = []
		yVals = []
		for i in range(len(concentrations)):
			y = list(combined.loc[elem, "Cal "+str(i+1)].as_matrix())
			x = [concentrations[i] for _ in range(len(y))]
			xVals += x
			yVals += y
		rWarnBool = True
		#Deals with cases where outliers exist
		slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(xVals), np.array(yVals))
		while rWarnBool:
			removed = False
			if np.abs(r_value) >= rWarnVal:
				rWarnBool = False
			elif np.abs(r_value) < rWarnVal:
				print("Warning: R for element " + elem + " is very low (" + str(r_value) + ")")
				print("Identifying outliers and reperforming fit...")
				#Separate rejected values
				stdev = np.std(yVals)
				mean = np.mean(yVals)
				holder_y = []
				holder_x = []
				rejects_y = []
				rejects_x = []
				j = 0
				for i in yVals:
					if i > mean-5*stdev and i < mean+5*stdev:
						holder_y.append(i)
						holder_x.append(xVals[j])
					else:
						rejects_y.append(i)
						rejects_x.append(xVals[j])
					j += 1
				xVals = holder_x
				yVals = holder_y
				print("Rejected Values: " + str(rejects_y) + str(rejects_x))
				removed = True
				slope, intercept, r_value, p_value, std_err = stats.linregress(np.array(xVals), np.array(yVals))
				#prevents re-iterating if problem not solved
				if removed:
					rWarnBool = False

			df = pd.DataFrame([slope, intercept, r_value, p_value, std_err], columns=[elem], index=["slope", "intercept", "r_value", "p_value", "std_err"])
		resultList += [df]
	return pd.concat(resultList, axis=1)

#gets the concentration values from a file if input, or manually from the user if necessary
def getConcentrations(num_cals, concFilePath, save_name="concentrations"):
    concentrations = None
    #if the user specified a file path, check if it's valid
    if concFilePath is not None:
        if os.path.exists(concFilePath):
            concFile = np.load(concFilePath)
            #make sure it has the right number of values
            if concFile.size != (num_cals):
                print ("input concentration file has wrong format. should have " + str(num_cals) + " rows")
            else:
                concentrations = concFile[:]
        else:
            print ("input concentration file not found")
    #if file input did't work, fallback to manual input
    if concentrations is None:
        concentrations = inputConcentrations(num_cals)
        np.save(save_name, concentrations)
    return concentrations

#find the dataframes for each .csv file
def Calibration_Script (file_input_path, output_dir):
    dataframeList, num_cals = getCalsFromDir(file_input_path)
    file_output_path = os.path.join(output_dir, "Calibration_Curves")

    if not os.path.exists(file_output_path):
        os.mkdir(file_output_path)

    if os.path.isfile(os.path.join(file_output_path,'concentrations.npy')):   	
    	bad_input = True
    	while bad_input:
    		use_cal = input('I see there is already a concentration file, would you like me to use it? [Y/N]: ')
    		if use_cal in ['Y', 'y']:
    			concFilePath = os.path.join(file_output_path,'concentrations.npy')
    			bad_input = False
    		elif use_cal in ['N', 'n']:
    			concFilePath = None
    			bad_input = False
    		else:
    			print("That's not an answer. Try again.")
    			use_cal = input('I see there is already a concentration file, would you like me to use it? [Y/N]: ')			
    else:
    	concFilePath = None

    if len(dataframeList) == 0:
        print ("Error: No .csv files found at directory " + file_input_path)
    else:
        #find the concentration values for each Cal
        concentrations = getConcentrations(num_cals, concFilePath, save_name=os.path.join(file_output_path, "concentrations"))

        #perform linear regression
        regResults = performLinearRegression(dataframeList, concentrations)

        #save linear regression results to a file
        regResults.to_csv(os.path.join(file_output_path, "CalibrationCurvesOut.csv"))
        #save cal columns to a csv file
        combined = pd.concat(dataframeList, axis=1)
        combined.to_csv(os.path.join(file_output_path, "CombinedCals.csv"))
        
    #return file paths for next step
    conc_file = os.path.join(file_output_path, "concentrations.npy")
    cal_file = os.path.join(file_output_path, "CalibrationCurvesOut.csv")
    cal_comb = os.path.join(file_output_path, "CombinedCals.csv")
    
    return cal_file, cal_comb, conc_file
