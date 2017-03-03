'''
Plotter.py
Written by J. Hayes; last editted 2/3/2017

Purpose: A script to plot the XRF spectra and their corresponding fits. Also plots the calibration curves

Notes: 
Uses pyplot from matplotlib; will ahave to include when building the executable
The main engine for the spectra plotting is Process_Fit_Spectra
The main engine for plotting the calibration curves is Process_Calibration_Curves
Main_Plotter runs both of the above in one shot (calibratin curves first, then XRF spectra)

Input:
1) calibration_dir - directory containing the calibration curve file and cal_concentratino file
2) cal_measure_dir - directory containing the combined calibration measurement files
3) spectra_dir - directory containing the XRF spectra and fit files
4) header_file - file containing the sample label information
5) output_dir - directory in which the output should be saved (folders inside this directory will be created)

Output:
Spectra and fits will be saved in output_dir/spectra_plots; plots for each plate will be saved in directories with the plate name
Calibration curves will be saved in calibration_dir/cal_curve_plots

'''

import pandas as pd 
import numpy as np
from os import listdir, mkdir
from os.path import isfile, join, isdir, basename
import matplotlib.pyplot as plt

#Make a directory for saving files in if it doesn't already exist; returns the full directory path
def Make_Save_Directory(output_dir, directory_name):
	directory = join(output_dir, directory_name)
	if not isdir(directory):
		mkdir(directory)
	return directory

#Relabel the columns in a spectra file to the sample names
def Relabel_Spectra_Columns(spectrum_file, plate_name, header_file):
	headers = pd.read_csv(header_file, index_col = 0)
	spectra = pd.read_csv(spectrum_file, index_col = False)
	plate_name = basename(spectrum_file).split('_')[0]

	for plate in headers.index:
		if plate.lower() == plate_name.lower():
			sample_names = headers.ix[plate, :]
			new_spectra_labels = []
			for name in sample_names:
				name = str(name)
				new_spectra_labels.append(name)
				new_spectra_labels.append(name)
				new_spectra_labels.append(name)
	
	#Try to relabel the spectra columns, if something goes wrong, the program will note it
	try:
		spectra.columns = new_spectra_labels
		spectra = Remove_Empty_Spectra(spectra)
	except UnboundLocalError:
		print("Hmm...something went wrong while labelling the spectra for plate %s" %plate_name)
		return
	
	return spectra

#Remove the spectra collected on the empty/blank cells from the spectra dataframe; takes dataframe, returns dataframe
def Remove_Empty_Spectra(data):
	empty_aliases = ['empty', 'Empty', 'x', 'X', 'blank', 'Blank']
	for column in set(data.columns): #Note: Use set because del data[column] will remove all instances of that column value	
		if column in empty_aliases:
			del data[column]
	return data

#Generate plots of each spectrum in a spectra file
def Spectra_Plotter(spectra, spectrum_file, output_dir):

	plate_name = spectrum_file.split('_')[0]
	save_dir = Make_Save_Directory(output_dir, plate_name)

	i = 0
	for column in spectra.columns:
		col_type = i % 3
		#Note range [2:] is result of first entry in column = energy, y_data, etc, second being blank
		if col_type == 0:
			energy = list(spectra.ix[:, i][2:])
		elif col_type == 1:
			y_data = list(spectra.ix[:, i][2:])
		elif col_type == 2:
			fit = list(spectra.ix[:, i][2:])
			Make_Spectrum_Plot(energy, y_data, fit, column, save_dir)
		i += 1

#Makes a single plot and saves it as output_dir/sample_name (note output_dir is whatever is specified when called)
def Make_Spectrum_Plot(energy, y_data, fit, sample_name, output_dir):
	file_name = join(output_dir, sample_name) + '.png'
	#Deal with sample repeat cases
	if isfile(file_name):
		i = 1 
		while isfile(file_name):
			file_name = file_name.replace('.png', '') + '_%s.png' %i
			i += 1 

	plt.plot(energy, y_data)
	plt.plot(energy, fit)
	plt.xlabel('Energy (keV)')
	plt.ylabel('Counts (a.u.)')
	plt.title(sample_name)
	print("Saving...%s" %file_name)
	plt.savefig(file_name)
	plt.clf()

#Script that processes the spectra files; master module for the spectra processing 
def Process_Fit_Spectra(spectra_dir, header_file, output_dir):
	spectra_files = [f for f in listdir(spectra_dir) if isfile(join(spectra_dir, f))]

	#Setup output directories if they don't already exist
	relabeled_spectra_dir = Make_Save_Directory(output_dir, "labeled_spectra")
	spectra_plot_dir = Make_Save_Directory(output_dir, "spectra_plots")

	for file in spectra_files:
		plate_name = file.split('_')[0]
		file_path = join(spectra_dir, file)

		print("Relabelling the plate: %s" %plate_name)
		spectra = Relabel_Spectra_Columns(file_path, plate_name, header_file)
		spectra.to_csv(join(relabeled_spectra_dir, file), index = False)
		
		print("Plotting spectra to %s" %spectra_plot_dir)		
		Spectra_Plotter(spectra, file, spectra_plot_dir)

#Read in the concentrations from the saved cal_concentrations.csv file
def Read_Concentrations(calibration_dir):
	conc_file = join(calibration_dir, "cal_concentrations.csv")
	conc_value_df = pd.read_csv(conc_file, header = None, index_col = 0)
	conc_values = np.array(conc_value_df[1])
	return conc_values

#Read in the calibration curve data, return at dictionary of tuples: Element: (slope, intercept)
def Read_Calibration_Curves(calibration_dir):
	calibration_curve_file = join(calibration_dir, "calibration_curves.csv")
	cal_curves = pd.read_csv(calibration_curve_file, index_col = 0)
	cal_curve_params = {}
	for element in cal_curves:
		cal_curve_params[element] = (cal_curves[element]['slope'], cal_curves[element]['intercept'])
	return(cal_curve_params)

#Read the averaged calibration measurements, return as dictionary or arrays: Element: [cal 1, cal 2, ...]
def Read_Calibration_Measurements(data_dir):
	cal_measurement_file = join(data_dir, 'calibration_data_avg.csv')
	cal_measurement_df = pd.read_csv(cal_measurement_file, index_col = 0)
	cal_measurements = {}
	for element in cal_measurement_df.index:
		cal_measurements[element] = np.array(cal_measurement_df.ix[element])
	return cal_measurements

#Plot the measured data and the fitted curve
def Make_Cal_Curve_Plot(cal_conc, cal_curve, cal_measurements, element, output_dir):
	file_name = join(output_dir, element + '.png')
	#print("Writing...%s" %file_name)
	#Calculate plot points for calibration curve
	
	cal_curve_y1 = cal_conc[0]*cal_curve[0] + cal_curve[1]
	cal_curve_y2 = cal_conc[-1]*cal_curve[0] + cal_curve[1]

	cal_measure_y = [cal_curve_y1, cal_curve_y2]
	cal_measure_x = [cal_conc[0], cal_conc[-1]]

	plt.plot(cal_measure_x, cal_measure_y)
	plt.scatter(sorted(cal_conc), sorted(cal_measurements))
	plt.xlabel('Concentration (ppm)')
	plt.ylabel('Peak Area (a.u.)')
	plt.title(element)
	plt.savefig(file_name)
	plt.clf()

#Batch process calibration data; main engine for plotting calibration curves and data
def Process_Calibration_Curves(calibration_dir, measurement_dir):
	output_dir = Make_Save_Directory(calibration_dir, "cal_curve_plots")
	calibration_concs = Read_Concentrations(calibration_dir)
	calibration_curves = Read_Calibration_Curves(calibration_dir)
	calibration_measurements = Read_Calibration_Measurements(measurement_dir)
	for element in calibration_measurements:
		Make_Cal_Curve_Plot(calibration_concs, calibration_curves[element], calibration_measurements[element], element, output_dir)

#The engine that runs entire script (spectra plotting and calibration curve plotting)
def Main_Plotter(calibration_dir, cal_measure_dir, spectra_dir, header_file, output_dir):
	Process_Calibration_Curves(calibration_dir, cal_measure_dir)
#	Process_Fit_Spectra(spectra_dir, header_file, output_dir)


header_file = 'C:/Users/John/Desktop/Test/labels.csv'
spectra_dir = 'C:/Users/John/Desktop/Test/Output/Spectra'
output_dir = 'C:/Users/John/Desktop/Test/Chibbar Group Work/XRF Analysis/XRF_Analysis/Non-Ca_Analysis/Old_Data/Output/'
calibration_dir = 'C:/Users/John/Desktop/Chibbar Group Work/XRF Analysis/XRF_Analysis/Non-Ca_Analysis/Old_Data/Output/calibration_results'
data_dir = 'C:/Users/John/Desktop/Chibbar Group Work/XRF Analysis/XRF_Analysis/Non-Ca_Analysis/Old_Data/Output/combined_files'

#Process_Fit_Spectra(spectra_dir, header_file, output_dir)
#Read_Concentrations(calibration_dir)
#df = Read_Calibration_Curves(calibration_dir)


Main_Plotter(calibration_dir, data_dir, spectra_dir, header_file, output_dir)
