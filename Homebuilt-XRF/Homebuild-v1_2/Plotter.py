import pandas as pd
import numpy as np
from os import listdir, mkdir
from os.path import isfile, join, isdir, split, dirname
import matplotlib.pyplot as plt

#Temporary for testing; the path where the spectra file to be plotted are saved
# spectra_dir = 'C:/Users/John/Desktop/Tester/Spectra'
# header_file = 'C:/Users/John/Desktop/Tester/Spectra/labels.csv'
# reg_info = 'C:/Users/John/Desktop/Tester/calibration/CalibrationCurvesOut.csv'
# cal_data = 'C:/Users/John/Desktop/Tester/calibration/CombinedCals.csv'
# cal_concs = 'C:/Users/John/Desktop/Tester/calibration/concentrations.npy'
# output_dir = 'C:/Users/John/Desktop/Tester/'
# calcium = False

#Rename the headers for all the files
def Rename_Spectra_Labels (spectra_dir, header_file):
	print("Working...one moment please")
	header = pd.read_csv(header_file, index_col = 0)
	all_files = [join(spectra_dir,f) for f in listdir(spectra_dir) if isfile(join(spectra_dir, f)) and 'ydata' in f]
	for f in all_files:
		if 'relabel' not in f and '.png' not in f:
			spec_data = pd.read_csv(f)
			for label in header.index.values:
				if label in f:
					i = 0
					for thing in header.loc[label]:
						spec_data = spec_data.rename(columns = {spec_data.columns[i]: str(thing)})
						spec_data = spec_data.rename(columns = {spec_data.columns[i+1]: str(thing) + ' data'})
						spec_data = spec_data.rename(columns = {spec_data.columns[i+2]: str(thing) + ' fit'})
						i += 3
					spec_data = spec_data.drop(spec_data.index[[0,1]])
			spec_data.to_csv(f.strip('.csv') + '_relabel.csv', index = False)

#Plot the spectra
def Plot_Spectra(spectra_dir, output_dir, calcium = False):
	if not isdir(join(output_dir, 'Spectra_Fits')):
		mkdir(join(output_dir, 'Spectra_Fits'))

	all_files = [f for f in listdir(spectra_dir) if isfile(join(spectra_dir, f)) and 'ydata' in f]
	for f in all_files:
		if 'relabel' in f:
			spec_data = pd.read_csv(join(spectra_dir,f))			
			for n in range(0,84, 3):
				x = list(spec_data.ix[:, n])
				y = list(spec_data.ix[:, n + 1])
				y_2 = list(spec_data.ix[:, n + 2])			
				plt.plot(x, y)
				plt.plot(x, y_2)
				#Change plotting window if it's calcium versus whole spectrum
				if calcium:
					plt.ylim(0,15000)
					plt.xlim(2.5,4.5)
				else:
					plt.xlim(2, 26)
				plt.xlabel('Energy (keV)')
				plt.ylabel('Counts (a.u.)')
				#Play around with this to better directory output rather than dumping all in one place...
				print("Now saving spectra fits")
				print("Saving..." + join(output_dir, 'Spectra_Fits', f.strip('relabel.csv')) + str(spec_data.columns[n]) + '.png')
				plt.savefig(join(output_dir, 'Spectra_Fits', f.strip('relabel.csv')) + str(spec_data.columns[n]) + '.png')
				plt.clf()

#Plot the calibration curves
def Plot_Calibrations(reg_info, cal_data, cal_concs, output_dir, calcium = False):
	reg_info = pd.read_csv(reg_info, index_col = 0)
	cal_data = pd.read_csv(cal_data, index_col = 0).transpose()
	concs = np.load(cal_concs) 

	if not isdir(join(output_dir, 'Calibration_Curves')):
		mkdir(join(output_dir, 'Calibration_Curves'))


	#Read in data and add calibration concentrations
	conc_list = []
	for sample in cal_data.index.values:
		for i in range(1,len(concs)+1):
			if ("Cal " + str(i)) in sample:
				conc_list.append(concs[i-1])
	cal_data['Concentration'] = conc_list

	#Make calibration curves then plot the curve with the data points 
	for element in reg_info.columns:
		if "element" not in element:
			x_vals = []
			y_vals = []
			if calcium:
				for i in range(0, 1101, 100):
					x_vals.append(i)
					y_vals.append(i*reg_info.ix['slope', element]+reg_info.ix['intercept', element])		
			else:
				for i in range(0, 111, 10):
					x_vals.append(i)
					y_vals.append(i*reg_info.ix['slope', element]+reg_info.ix['intercept', element])

			plt.scatter(cal_data['Concentration'], cal_data[element], color = 'k')
			plt.plot(x_vals, y_vals, ls = '--', color = '0.75', lw = 2)
			
			#Set the plotting window
			if calcium:
				plt.xlim(0,1100)
			else:
				plt.xlim(0,110)
			
			#Window dressing
			plt.xlabel('Concentration (ppm)')
			plt.ylabel('Normalized Area (a.u.)')

			#Save the figure; probably want to change to include output directory
			print("Now saving calibration curves")
			print("Saving..." + join(output_dir, 'Calibration_Curves', str(element)) + '.png')
			plt.savefig(join(output_dir, 'Calibration_Curves', str(element)) + '.png')
			plt.clf()


def Plotter_Program(spectra_dir, header_file, output_dir, cal_data, calcium, cal_concs, reg_info):
	Rename_Spectra_Labels(spectra_dir, header_file)
	Plot_Spectra(spectra_dir, output_dir, calcium)
	Plot_Calibrations(reg_info, cal_data, cal_concs, output_dir, calcium)