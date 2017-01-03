import FormatPlateData
import ConcentrationCalc
import CalibrationCurves
import NormalizeData
import Plotter
from os.path import isfile, join, exists, abspath, isdir
from os import mkdir, remove

def dir_reader():
	with open('dir_parameters.txt', 'r') as f:
		lines = f.readlines()
		fit_dir = lines[0].strip('\n')
		header_file = lines[1].strip('\n')
		io_dir = lines[2].strip('\n')
		output_dir = lines[3].strip('\n')
		is_cal = lines[4].strip('\n')
		
		if is_cal == "False":
			is_cal = False
		elif is_cal == "True":
			is_cal = True

	print("\n Here are the files to input: ")
	print("Fit Directory: " + fit_dir)
	print("Header file: " + header_file)
	print("Io Directory: " + io_dir)
	print("Output Directory " + output_dir)
	print("Calcium Analysis: " + str(is_cal))
	

	is_correct = input('Do you want to keep these parameters? [Y/N]')
	valid_response = False
	while not valid_response:
		if is_correct in ['Y', 'y']:
			valid_response = True
			return fit_dir, header_file, io_dir, output_dir, is_cal
		elif is_correct in ['N', 'n']:
			fit_dir, header_file, io_dir, output_dir, is_cal = File_Directorys()
			return fit_dir, header_file, io_dir, output_dir, is_cal
		else:
			print("Sorry, I don't understand...")
			is_correct = input('Do you want to keep these parameters? [Y/N]')

def dir_parameters_exists():
	if isfile('dir_parameters.txt'):
		bad_input = True
		while bad_input:
			user_param = input('I noticed I happened to write down the previous directory parameters. Would you like me to use them? [Y/N]:  ')
			if user_param in ['Y', 'y']:
				dir_parameters = True
				bad_input = False
				return dir_parameters
			elif user_param in ['N', 'n']:
				print("That's okay, I'll ask for new values.")
				bad_input = False
			else:
				print("I don't understand. Please try again.")
				user_param = input("Would you like me to use the previous directory parameters? [Y/N]")

def Get_Fit_Dir():
	fit_dir = input("Input folder containing the fits to be analyzed: ")
	
	while not isdir(fit_dir):
		print("Not a valid path.")
		fit_dir = input("Input folder containing the fits to be analyzed: ")

	return fit_dir

def Get_Output_Dir():
	output_dir = input("Where do you want to save the analysis: ")
	
	while not isdir(output_dir):
		create_dir = input(output_dir + " does not exist, would you like me to create it? [Y/N]   : ")
		if create_dir in ['Y','y']:
			mkdir(output_dir)
		elif create_dir in ['N','n']:
			output_dir = input("Where do you want to save the analysis: ")
		else:
			print("I'm sorry, I don't understand")

	return output_dir

def Get_Io_Dir():
	io_dir = input("Input folder containing the Io files: ")
	
	while not isdir(io_dir):
		print("Not a valid path.")
		io_dir = input("Input folder containing the Io files: ")

	return io_dir

def Get_Header_File():
	header_file = input("Header file (with path): ")
	
	while not isfile(header_file):
		print("Not a valid file.")
		header_file = input("Header file (with path): ")

	return header_file

def Ca_Analysis_Question():
	calcium = input('Is this the Calcium Analysis? [Y/N]: ')
	valid_response = False
	while not valid_response: 
		if calcium in ['Y' , 'y']:
			calcium = True
			valid_response = True
			return calcium
		elif calcium in ['N' , 'n'] or calcium == None:
			calcium = False
			valid_response = True
			return calcium 
		else:
			calcium = input("Sorry, I didn't get that. Is this the Calcium Analysis? (Y/N; Default is N): ")

def Action_Decider(fit_dir, header_file, io_dir, output_dir, calcium):
	whole_analysis = False
	while not whole_analysis:
		print("What would you like to do? Type in the number to indicate procedure: ")
		print("1) Run full analysis")
		print("2) Format plate data")
		print("3) Normalize data and add headers")
		print("4) Perform Linear Regression")
		print("5) Calculate concentrations and perform error analysis")
		print("6) Generate plots of the fits and calibration curves")
		print("7) Exit the program without deleting directory parameters file")
		print("8) Exit program and delete directory parameters file")
		action = input('Your choice: ')
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
