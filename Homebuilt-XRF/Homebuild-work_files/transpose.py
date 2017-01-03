import pandas as pd 
import numpy as np 
from os import listdir, walk, mkdir
from os.path import isfile, join, exists, abspath, isdir, dirname, realpath
import sys

def DataFinder():
	path = input("Please input the data directory path: ")
	
	while not isdir(path):
		print("Sorry, I didn't get that...")
		path = input("Please input the data directory path: ")

	return path

def Transpose(datafiledir):
	if not exists(join(datafiledir, 'Transposed')):
		mkdir(join(datafiledir, 'Transposed'))
	
	out_dir = join(datafiledir, 'Transposed')
	
	all_files = [f for f in listdir(datafiledir) if isfile(join(datafiledir, f))]
	for f in all_files:
		print(join(out_dir, f.strip('.dat')) + '_trans.csv')
		if 'Ge13El' in f and 'Raw' not in f and 'Hline' not in f and 'Vline' not in f and 'trans' not in f:
			data = pd.read_csv(join(datafiledir, f), delimiter = '\t', header = None, index_col = False)
			data = data.transpose()
			data.to_csv(join(out_dir, f.strip('.dat')) + '_trans.csv', index = False, header = False)

datafiledir = DataFinder()
Transpose(datafiledir)

	

