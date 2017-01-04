import pandas as pd 
from collections import Counter
import sys

#Reads a dataframe and returns a parsed dataframe with only the specified
#taxonomy keys (a way to break down the big list into smaller lists by family)
def Taxonomy_Parse(df, tax_key):
	data_dict = {}
	for row in df.iterrows():
		taxonomy = row[1]["Taxonomy Info"]
		tax_list = taxonomy.split(',')
		
		#Remove gratitous spaces from list
		new_tax_list = []
		for item in tax_list:
			item = item.strip(' ')
			new_tax_list.append(item)
		tax_list = new_tax_list
		
		#Parse the taxonomy list, if a hit, add to dictionary
		if tax_key in tax_list:
			key = row[1]['Accession Number']
			species = row[1]['Species']
			description = Description_Stripper(row[1]['Description/Name'])
			data_dict[key] = [species, description, taxonomy]

	if len(data_dict) != 0:
		df2 = pd.DataFrame.from_dict(data_dict, orient = 'index')
	else:
		print("Whoops, something went wrong")
		df2 = pd.DataFrame()

	return df2 

#Removes species information from the description section
def Description_Stripper(description):
	if '[' in description:
		record = True
		new_name = ''
		for char in description:
			if char == '[':
				record = False
			if record:
				new_name += char
	else:
		new_name = description
	return new_name


if len(sys.argv) > 2:
	for tax_key in sys.argv[2:]:
		df = pd.read_csv(sys.argv[1])
		Taxonomy_Parse(df, tax_key).to_csv('%s.csv' %tax_key, header = False)

