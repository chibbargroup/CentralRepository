from os.path import isfile, join, isdir
from os import mkdir, remove
from User_Input import Get_Input_Parameters, Yes_No_Response
from FormatPlateData import ExtractDirectory
from NormalizeData import Batch_Process


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
	valid_choices = ['1', '2', '3', '4', '5', '6', '7'] #Note: This needs to be updated as features/options are added
	valid_response = False
	while not valid_response:
		action = input('Your choice: ')
		if action in valid_choices:
			valid_response = True
		else:
			print("Sorry, that doesn't seem to work for me. Try again?")
	return action

#Runs the FormatPlateData script; doesn't need to be a module, but here for consistency
def Format_Plate_Data(fit_dir, output_dir):
	ExtractDirectory(fit_dir, output_dir)

def Normalize_Data(header_file, io_file_dir, output_dir):
	ca_analysis = Yes_No_Response("Is this a Ca analysis? [Y/N]: ")
	peakfit_dir = join(output_dir, 'peakfits')
	Batch_Process(peakfit_dir, header_file, io_file_dir, output_dir, ca_analysis)

process_monitor = {'FormatPlate': False, 'Normalize': False, 'FileCombine': False,
                   'LinRegress': False, 'ConcCalc': False}

print(process_monitor)

#Prints out possible actions the user can take; needs to be updated as features are added/removed 
def Action_Decider():
	print("What would you like to do? Type in the number to indicate procedure: ")
	print("1) Run full analysis")
	print("2) Format plate data")
	print("3) Normalize data and add headers")
	print("4) Perform Linear Regression")
	print("5) Calculate concentrations and perform error analysis")
	print("6) Generate plots of the fits and calibration curves")
	print("7) Exit the program without deleting directory parameters file")
	print("8) Exit program and delete directory parameters file")
	action = Valid_Decision()
	return action






'''



def Action_Decider(fit_dir, header_file, io_dir, output_dir, calcium):
	whole_analysis = False
	while not whole_analysis:

		action = int(action)
		if action == 1:
			FormatPlateData.extractDirectory(fit_dir, output_dir)
			NormalizeData.Replace_Header_Batch(join(output_dir,'peakfits'), io_dir, header_file, output_dir, calcium)
			cal_file, cal_comb, conc_file = CalibrationCurves.Calibration_Script(join(output_dir,'Norm'), output_dir)
			ConcentrationCalc.Analyzer(join(output_dir,'Norm'), cal_file, cal_comb, conc_file, output_dir)
		elif action == 2:
			FormatPlateData.extractDirectory(fit_dir, output_dir)
		elif action == 3:
			NormalizeData.Replace_Header_Batch(join(output_dir,'peakfits'), io_dir, header_file, output_dir, calcium)
		elif action == 4:
			cal_file, cal_comb, conc_file = CalibrationCurves.Calibration_Script(join(output_dir,'Norm'), output_dir)
		elif action == 5:
			ConcentrationCalc.Analyzer(join(output_dir,'Norm'), cal_file, cal_comb, conc_file, output_dir)
		elif action == 6:
			Plotter.Plotter_Program(join(output_dir, 'spectra'), header_file, output_dir, cal_comb, calcium, conc_file, cal_file)
		elif action == 7:
			whole_analysis = True
			dir_parameters_writer(fit_dir, header_file, io_dir, output_dir, calcium)
			break
		elif action == 8:
			remove('dir_parameters.txt')
			whole_analysis = True
			break
		else:
			action = input("I'm sorry Dave, I can't let you do that. Please input a number between 1-6: ")

def File_Directorys():
	is_right = False
	while not is_right:
		fit_dir = Get_Fit_Dir()
		header_file = Get_Header_File()
		io_dir = Get_Io_Dir()
		output_dir = Get_Output_Dir()
		calcium = Ca_Analysis_Question()

		print("\n Here are the files to input: ")
		print("Fit Directory: " + fit_dir)
		print("Header file: " + header_file)
		print("Io Directory: " + io_dir)
		print("Output Directory " + output_dir)
		print("Calcium Analysis: " + str(calcium))

		right = input("Is this correct [Y/N]: ")
		if right in ['Y', 'y']:
			is_right = True
			dir_parameters_writer(fit_dir, header_file, io_dir, output_dir, calcium)
		elif right in ['N', 'n']:
			is_right = False
		else:
			print("I don't get it...can you try again?")
			right = input("Is this correct [Y/N]: ")

	return fit_dir, header_file, io_dir, output_dir, calcium

def dir_parameters_writer(fit_dir, header_file, io_dir, output_dir, is_calcium):
	with open('dir_parameters.txt', 'w') as f:
		f.write(fit_dir + '\n')
		f.write(header_file + '\n')
		f.write(io_dir + '\n')
		f.write(output_dir + '\n')
		f.write(str(is_calcium))


dir_parameter = dir_parameters_exists()
if dir_parameter:
	fit_dir, header_file, io_dir, output_dir, is_calcium = dir_reader()
else:
	fit_dir, header_file, io_dir, output_dir, is_calcium = File_Directorys()
	dir_parameters_writer(fit_dir, header_file, io_dir, output_dir, is_calcium)
Action_Decider(fit_dir, header_file, io_dir, output_dir, is_calcium)
'''