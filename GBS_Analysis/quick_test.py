from collections import OrderedDict

file = './loc_test.tab'
test_dict = OrderedDict()

with open(file, 'r') as f:
	f.readline()
	for line in f:
		split_line = line.split('\t')
		key = split_line[0]
		data = '\t'.join(split_line[1:])
		test_dict[line] = data

print(list(test_dict.values()))