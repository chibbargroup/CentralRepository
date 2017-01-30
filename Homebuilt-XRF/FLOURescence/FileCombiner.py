'''
FileCombiner.py
Written by J. Hayes, Last edited 1/27/2017

Purpose: Combines the normalized data into a single data frame and then parses out
the calibration data and sample data into separate files. Then determines the average
value for each calibration point and sample point

Inputs:
1) norm_data_file_dir: Directory which contains the normalized peak area data
2) output_dir: Where the resulting files should be saved

Output:
Program saves 4 files:
1) calibration_data.csv: the raw, unaveraged calibration data
2) calibration_data_avg.csv: the averaged calibration data
3) sample_data.csv: the raw, unaveraged sample data
4) sample_data_avg.csv: the averaged sample data

'''

import pandas as pd
import numpy as np
import collections
from os import listdir
from os.path import isfile, join

#Define OrderedSet class for data averaging
class OrderedSet(collections.MutableSet):
    def __init__(self, iterable=None):
        self.end = end = [] 
        end += [None, end, end]         # sentinel node for doubly linked list
        self.map = {}                   # key --> [key, prev, next]
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            end = self.end
            curr = end[1]
            curr[2] = end[1] = self.map[key] = [key, curr, end]

    def discard(self, key):
        if key in self.map:        
            key, prev, next = self.map.pop(key)
            prev[2] = next
            next[1] = prev

    def __iter__(self):
        end = self.end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = self.end[1][0] if last else self.end[2][0]
        self.discard(key)
        return key

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return set(self) == set(other)

#Open all data files and combine into single dataframe
def Data_Combiner(norm_file_dir):
    data = pd.DataFrame()
    all_files = all_files = [f for f in listdir(norm_file_dir) if isfile(join(norm_file_dir, f))]
    for f in all_files:
        df = pd.read_csv(join(datafiledir,f), index_col = 0)
        data = pd.concat([data,df], axis = 1)
    data = data.drop('Io')
    return data 

#Parse through combined data and pick out data columns
def Sample_Data_Parse(data):
	empty_aliases = ['empty', 'Empty', 'x', 'X', 'blank', 'Blank'] #Can be modified as necessary
	for column in set(data.columns): #Note: Use set because del data[column] will remove all instances of that column value
		column_read = column.split('.')[0]	
		if column_read in empty_aliases:
			del data[column]
		elif "Cal " in column:
			del data[column]
	return data

#Parse through combined data nd pick out calibration sample columns
def Cal_Data_Parse(data):
	for column in set(data.columns):
		if "Cal " not in column:
			del data[column]
	return data

#Rename columns to have format sample.#
def Replicate_Renamer(data):
    column_list = list(data.columns)
    new_names = []
    for column in column_list:
        if '.' not in column:
            column = column + '.1'
        if column not in new_names:
            new_names.append(column)
        else:
            i = 1
            rename = True
            while rename:
                if '.' in column:
                    column = column.split('.')[0]               
                column = column + '.%s' %i
                if column not in new_names:
                    new_names.append(column)
                    rename = False
                else:
                    i += 1
    data.columns = new_names
    return data

#Average the data
def Data_Averager(data):
    averaged_data = pd.DataFrame()
    
    #Make list of each individual sample
    sample_list = []
    for column in data:
        column = column.split('.')[0]
        sample_list.append(column)
    sample_list = list(OrderedSet(sample_list))

    #Average all the data columns containing data from that sample 
    for sample in sample_list:
        column_positions = [i for i,j in enumerate(data.columns) if sample == j.split('.')[0]]
        sample_data = data.iloc[:,column_positions]
        sample_average = sample_data.sum(axis=1)/len(column_positions)
        sample_average = sample_average.rename(sample)
        averaged_data = pd.concat([averaged_data, sample_average], axis = 1)
    return(averaged_data)

#Combined and average sample data
def Sample_Data_Separator(norm_data_dir, output_dir):
	df = Data_Combiner(norm_data_dir)
	norm_data = Sample_Data_Parse(df)
	norm_data = Replicate_Renamer(norm_data)
	norm_data_avg = Data_Averager(norm_data)
	norm_data.to_csv(join(output_dir, 'sample_data.csv'))
	norm_data_avg.to_csv(join(output_dir, 'sample_data_avg.csv'))

#Combine and average calibration data
def Calibration_Data_Separator(norm_data_dir, output_dir):
	df = Data_Combiner(norm_data_dir)
	cal_samples = Cal_Data_Parse(df)
	cal_samples = Replicate_Renamer(cal_samples)
	cal_samples_avg = Data_Averager(cal_samples)
	cal_samples.to_csv(join(output_dir, 'calibration_data.csv'))
	cal_samples_avg.to_csv(join(output_dir, 'calibration_data_avg.csv'))

#Run the Script
def Data_Separator(norm_data_dir, output_dir):
	Sample_Data_Separator(norm_data_dir, output_dir)
	Calibration_Data_Separator(norm_data_dir, output_dir)

