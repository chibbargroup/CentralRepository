from os.path import isfile, join, isdir
from os import mkdir, remove

#Asks a Yes/No question and returns either True (yes) or False (no) based on the user response
def Yes_No_Response(question_string):
	positive_response = ['yes', 'y', 'yup', 'yea', 'yeah']
	negative_response = ['no', 'n', 'nope', 'nah', 'nay']
	valid_response = False
	while not valid_response:
		response = input(question_string)
		if response.lower() in positive_response:
			valid_response = True
			return True
		elif response.lower() in negative_response:
			valid_response = True
			return False
		else:
			print("I'm sorry Dave, I can't do that. Please input a valid response.")

#Get the file containing the header information; file must have the .csv file extension
def Header_File_Input():
	valid_file = False
	while not valid_file:
		header_file = input("Please direct me to the .csv file containing the sample label information: ")
		print(header_file[-4:])
		if isfile(header_file) and header_file[-4:] == '.csv':
			valid_file = True
		else:
			print("You seem to have me HEADED in the wrong direction. I can't find that file...")
	return header_file

#Get the directory containing the FIT files/folders
def Fit_Directory_Input():
	valid_dir = False
	while not valid_dir:
		fit_directory = input("Please direct me to the folder containing the PyMCA Fit files: ")
		if isdir(fit_directory):
			valid_dir = True
		else:
			print("Hmm...something doesn't FIT. I can't find that directory...")
	return fit_directory

#Direct the script to an output directory; if it doesn't exist, then see if the user would like to create it
def Output_Directory_Input():
	valid_dir = False
	while not valid_dir:
		fit_directory = input("Where would you like me to save the results? Please input a location: ")
		if isdir(fit_directory):
			valid_dir = True
		else:
			make_new_dir = Yes_No_Response("That directory doesn't seem to exist, would you like me to create it? [Y/N]: ")
			if make_new_dir:
				mkdir(fit_directory)
				valid_dir = True
			else:
				print("I'm not sure what you want to get OUT of this relationship...")
	return fit_directory

#Gets the fit directory, header file, and output directory via user input 
def User_Input():
	fit_dir = Fit_Directory_Input()
	header_file = Header_File_Input()
	output_dir = Output_Directory_Input()
	return fit_dir, header_file, output_dir

#Writes the parameters used for the analysis to the file prev_params.txt file in the script directory
def Previous_Parameters_Writer(fit_dir, header_file, output_dir):
	with open('prev_params.txt', 'w') as f:
		f.write(fit_dir + '\n')
		f.write(header_file + '\n')
		f.write(output_dir + '\n')

#Reads the file prev_params.txt, which is by default stored in the script directory
def Previous_Parameter_Reader():
	parameters = []
	with open('prev_params.txt', 'r') as f:
		for line in f:
			parameters.append(line.strip('\n'))
	fit_dir = parameters[0]
	header_file = parameters[1]
	output_dir = parameters[2]
	return fit_dir, header_file, output_dir

#Gets the fit directory, header file, and output directory, either by reading prev_params.txt, or User_Input
def Parameter_Getter():
	if isfile('prev_params.txt'):
		use_previous = Yes_No_Response("I see there are previous input parameters, would you like me to use them? ")
		if use_previous:
			fit_dir, header_file, output_dir = Previous_Parameter_Reader()
		else:	
			fit_dir, header_file, output_dir = User_Input()
	else:
		fit_dir, header_file, output_dir = User_Input()

	return fit_dir, header_file, output_dir

#Runs Parameter_Getter and runs until the user is happy with the parameters, if yes, will parameters to prev_params.txt
def Get_Input_Parameters():
	fit_dir, header_file, output_dir = Parameter_Getter()
	user_satisfied = False
	while not user_satisfied:
		print("Okay, here's what I have for the analysis: ")
		print("The specified fit directory is: %s" %fit_dir)
		print("The specified header/sample label file is: %s" %header_file)
		print("The specified output directory is: %s" %output_dir)
		
		correct_parameters = Yes_No_Response("Is this correct? ")
		if correct_parameters:
			user_satisfied = True
		else:
			print("Alright, that's cool. Let's try this again then")
			fit_dir, header_file, output_dir = User_Input()
	Previous_Parameters_Writer(fit_dir, header_file, output_dir)
	return fit_dir, header_file, output_dir