#Local imports related to FLOURescence only
from UserInput import Get_Input_Parameters, Yes_No_Response
from FormatPlateData import ExtractDirectory
from NormalizeData import Batch_Process
from FileCombiner import Data_Separator
from LinearRegression import Linear_Regression_Script
from ConcentrationCalculator import Concentration_Analysis
from PlotFigures import Main_Plotter
#Global Library Imports
from os.path import isfile, join, isdir
from os import mkdir, remove
import sys


#Print the program boiler plate/header; to be run upon initializing the script
#Make this more elegant in future revs
def Boiler_Plate_Printer():
	print("***********************************************")
	print("*                 Welcome to                  *")
	print("*              FLOURescence v1.0              *")
	print("*          By D. Sanche and J. Hayes          *")
	print("*                (C) Jan 2017                 *")
	print("***********************************************")

#Takes user input, determines if it is valid/can be read or not; if not, forces user to give good input
def Valid_Decision():
	valid_choices = ['1', '2', '3', '4', '5', '6', '7', '8'] #Note: This needs to be updated as features/options are added
	valid_response = False
	while not valid_response:
		action = input('Your choice: ')
		if action in valid_choices:
			valid_response = True
		else:
			print("Sorry, that doesn't seem to work for me. Try again?")
	return action

#Runs the FormatPlateData script; doesn't need to be a module, but here for consistency
def Format_Plate_Data(fit_dir, output_dir, process_monitor):
	ExtractDirectory(fit_dir, output_dir)
	process_monitor['FormatPlate'] = True
	return process_monitor

#Runs the normalize data script; saves the output in output_dir/normalized
def Normalize_Data(header_file, io_file_dir, output_dir, process_monitor):
	ca_analysis = Yes_No_Response("Is this a Ca analysis? [Y/N]: ")
	peakfit_dir = join(output_dir, 'peakfits')
	Batch_Process(peakfit_dir, header_file, io_file_dir, output_dir, ca_analysis)
	process_monitor['Normalize'] = True
	return process_monitor

#Runs the file combiner script, saves the output in output_dir/combined_files
def File_Combiner(output_dir, process_monitor):
	norm_file_dir = join(output_dir, 'normalized')
	new_output_dir = join(output_dir, 'combined_files')
	if not isdir(new_output_dir):
		mkdir(new_output_dir)
	Data_Separator(norm_file_dir, new_output_dir)
	process_monitor['FileCombine'] = True
	return process_monitor

#Runs the linear regression script, saves the output in output_dir/calibration_results
def Linear_Regression(output_dir,process_monitor):
	data_dir = join(output_dir, 'combined_files')
	cal_data_file = join(data_dir, 'calibration_data_avg.csv')
	new_output_dir = join(output_dir, 'calibration_results')
	if not isdir(new_output_dir):
		mkdir(new_output_dir)
	Linear_Regression_Script(cal_data_file, new_output_dir)
	process_monitor['LinRegress'] = True
	return process_monitor

#Runs the concentration calculator script; saves the output in output_dir
def Concentration_Calc(output_dir, process_monitor):
	data_dir = join(output_dir, 'combined_files')
	calibration_dir = join(output_dir, 'calibration_results')
	raw_sample_data = join(data_dir, 'sample_data.csv')
	avg_sample_data = join(data_dir, 'sample_data_avg.csv')
	cal_curves = join(calibration_dir, 'calibration_curves.csv')
	Concentration_Analysis(raw_sample_data, avg_sample_data, cal_curves, output_dir)
	process_monitor['ConcCalc'] = True
	return process_monitor

#Generates plots of the spectra and their fits as well as the calibration curves; not yet implemented
def Generate_Plots(output_dir, header_file, process_monitor):
	calibration_dir = join(output_dir, 'calibration_results')
	cal_measurement_dir = join(output_dir, 'combined_files')
	spectra_dir = join(output_dir, 'spectra')
	Main_Plotter(calibration_dir, cal_measurement_dir, spectra_dir, header_file, output_dir)
	process_monitor['Plots'] = True
	return process_monitor

#Prints out possible actions the user can take; needs to be updated as features are added/removed 
def Action_Decider():
	print("What would you like to do? Type in the number to indicate procedure: ")
	print("1) Format plate data")
	print("2) Normalize the data and add sample labels")
	print("3) Combine the normalized data")
	print("4) Perform the linear regresssion analysis")
	print("5) Calculate concentrations and perform error analysis")
	print("6) Generate plots of the fits and calibration curves") #Not yet implemented
	print("7) Do all of the above (this may overwrite existing data in the output directory)")
	print("8) Exit the program")
	action = Valid_Decision()
	return action

#Reads the action decided by the user and runs it
#Note: Would be cleaner to run the try anyway functionality as a separate function, but need to implement/figure out how to do that first
def Action_Taker(action, fit_dir, header_file, io_file_dir, output_dir, process_monitor):
	action = int(action)
	if action == 1:
		process_monitor = Format_Plate_Data(fit_dir, output_dir, process_monitor)
		return process_monitor

	if action == 2 and process_monitor['FormatPlate']:
		process_monitor = Normalize_Data(header_file, io_file_dir, output_dir, process_monitor)
		return process_monitor
	elif action == 2:
		print("It seems you haven't run step 1 yet. Running this step right now may break me.")
		try_anyway = Yes_No_Response("Would you like to try to run this anyways? [Y/N] ")
		if try_anyway:
			try:
				process_monitor = Normalize_Data(header_file, io_file_dir, output_dir, process_monitor)
				return process_monitor
			except (ValueError, FileNotFoundError, OSError):
				print("Yup, you broke me :( try running the previous step first")
				return process_monitor

	if action == 3 and process_monitor['Normalize']:
		process_monitor = File_Combiner(output_dir, process_monitor)
		return process_monitor
	elif action == 3:
		print("It seems you haven't run step 2 yet. Running this step right now may break me.")
		try_anyway = Yes_No_Response("Would you like to try to run this anyways? [Y/N] ")
		if try_anyway:
			try:
				process_monitor = File_Combiner(output_dir, process_monitor)
				return process_monitor
			except (ValueError, FileNotFoundError, OSError):
				print("Yup, you broke me :(...try running the previous step first")
				return process_monitor

	if action == 4 and process_monitor['FileCombine']:
		process_monitor = Linear_Regression(output_dir, process_monitor)
		return process_monitor
	elif action == 4:
		print("It seems you haven't run step 3 yet. Running this step right now may break me.")
		try_anyway = Yes_No_Response("Would you like to try to run this anyways? [Y/N] ")
		if try_anyway:
			try:
				process_monitor = Linear_Regression(output_dir, process_monitor)
				return process_monitor
			except (ValueError, FileNotFoundError, OSError):
				print("Yup, you broke me :(...try running the previous step first")
				return process_monitor

	if action == 5 and process_monitor['LinRegress']:
		process_monitor = Concentration_Calc(output_dir, process_monitor)
		return process_monitor
	elif action == 5:
		print("It seems you haven't run step 4 yet. Running this step right now may break me.")
		try_anyway = Yes_No_Response("Would you like to try to run this anyways? [Y/N] ")
		if try_anyway:
			try:
				process_monitor = Concentration_Calc(output_dir, process_monitor)
				return process_monitor
			except (ValueError, FileNotFoundError, OSError):
				print("Yup, you broke me :(...try running the previous step first")
				return process_monitor

	if action == 6 and process_monitor['LinRegress'] and process_monitor['FileCombine']:
		process_monitor = Generate_Plots(output_dir, header_file, process_monitor)
		return process_monitor
	elif action == 6:
		print("It seems you haven't run the required previous steps yet. Running this step right now may break me.")
		try_anyway = Yes_No_Response("Would you like to try to run this anyways? [Y/N] ")
		if try_anyway:
			try:
				process_monitor = Generate_Plots(output_dir, header_file, process_monitor)
				return process_monitor
			except (ValueError, FileNotFoundError, OSError):
				print("Yup, you broke me :(...try running the previous step first")
				return process_monitor
				
	if action == 7:
		Format_Plate_Data(fit_dir, output_dir, process_monitor)
		Normalize_Data(header_file, io_file_dir, output_dir, process_monitor)
		File_Combiner(output_dir, process_monitor)
		Linear_Regression(output_dir, process_monitor)
		Concentration_Calc(output_dir, process_monitor)
		Generate_Plots(output_dir, header_file, process_monitor)
		return process_monitor

	elif action == 8:
		sys.exit()
	return process_monitor

def Show_Runner():
	Boiler_Plate_Printer()
	fit_dir, header_file, io_file_dir, output_dir = Get_Input_Parameters()
	process_monitor = {'FormatPlate': False, 'Normalize': False, 'FileCombine': False,
					'LinRegress': False, 'ConcCalc': False, 'Plots': False}
	stop = False
	while not stop:
		action = Action_Decider()
		process_monitor = Action_Taker(action, fit_dir, header_file, io_file_dir, output_dir, process_monitor)

Show_Runner()

