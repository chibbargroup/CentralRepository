#OneMap_File_Generator_outcross.py
#This script will generate from a collection of vcf SNP calls a input OneMap_File_Generator
#for OneMap in order to facilitate linkage mapping. For outcross analysis only.
#Input: A file containing a list of file locations of the files to be analyized.
#Output: Onemap Input file in the proper format.



#main method for the script
def main(inputVCFFile,inputParent1Loc = '',inputParent2Loc = ''):
  import os.path
  
  #check that file exits
  if os.path.isfile(inputVCFFile) == False:
    print('Specified File does not exist: ' + inputVCFFile)
    return
  
  
  
  #initialize data var
  markerDict = {}
  statsDict = {}
  fileList = list()   
  knownParent = True
   
  LoadVCFTabFile(inputVCFFile,markerDict,statsDict,fileList)    
    
  #if inputParent1Loc and inputParent2Loc are not set run a unknown parent analysis
  if inputParent1Loc == '' and inputParent2Loc == '':
    knownParent = False
    p1_index = -1
    p2_index = -1
  else:    
    knownParent = True
    p1_index = fileList.index(inputParent1Loc)
    p2_index = fileList.index(inputParent2Loc)  
    
  numOfMarkers = 0  
  outFileString = ''
  
  
    
  for key in markerDict:
    if knownParent:
      markerTypeData = DeterminMarkerType_knownParent(p1_index,p2_index,markerDict[key])
    else: #unknown parent
      markerTypeData = DeterminMarkerType_unknownParent(markerDict[key])

      
    if markerTypeData['Type'] == 'F': #failed marker determination
      continue
    markerEntry = GenerateMarkerEntry(key,markerTypeData,markerDict[key],p1_index,p2_index)
    if markerEntry != '':
      numOfMarkers += 1
      outFileString += markerEntry + '\n\n'
      
      
  if knownParent:
    numOfSamples = len(fileList)-2   
  else:
    numOfSamples = len(fileList)   
  
  
  
  #write out file
  outFile = open(inputVCFFile + '.onemap','w')
  
  outFile.write(str(numOfSamples) + ' ' + str(numOfMarkers) + '\n')
  outFile.write(outFileString)
  
  outFile.close()
  
  
   
   
   
#Load a tab version of the VCF file and create a dictionary
#of markers and return it to calling function.
def LoadVCFTabFile(filePath,markerDict,statsDict,fileList):
  
  inputFile = open(filePath,'r')
  
  
  for line in inputFile:
    if '#CHROM' in line: #header line, fill out file list
      split_line = line.split('\t')
      for section in split_line:
	if section == '#CHROM' or section == 'POS' or section == 'REF':
	  continue
	else:
	  fileList.append(section)
    else:
      split_line = line.split('\t')
      
      #key name will be #Chrom + _p + pos
      keyName = split_line[0] + '_p' + split_line[1];
      
      #StatsDict format: snp position, # of missing data
      statsDict[keyName] = list()
      statsDict[keyName].append(split_line[1])
      statsDict[keyName].append(0)
      
      markerDict[keyName] = list()
      
      for i in range(3,len(split_line)):
	markerDict[keyName].append(split_line[i].strip('\n'))
	if split_line[i] == './.': #missing data
	  statsDict[keyName][1] += 1
	  
  inputFile.close()
	  

#determins the marker type based on Linkage paper
#Rongling Wu, Chang-Xing Ma, Ian Painter, Zhao-Bang Zeng, Simultaneous Maximum Likelihood Estimation of Linkage and Linkage Phases in #Outcrossing Species, Theoretical Population Biology, Volume 61, Issue 3, May 2002, Pages 349-363, ISSN 0040-5809, #http://dx.doi.org/10.1006/tpbi.2002.1577.
#(http://www.sciencedirect.com/science/article/pii/S0040580902915777)
#IF THEN ELSE chain in function derived from table 1 in the refrenced paper.  There's got to be a better way to figure this out.
def DeterminMarkerType_knownParent(p1_idx,p2_idx,allMarkerData):
  
  selectedType = ''
  a = ''
  b = ''
  c = ''
  d = ''
  
  p1Marker_s = allMarkerData[p1_idx].split('/')
  p2Marker_s = allMarkerData[p2_idx].split('/')
  
  
  p1_a1 = p1Marker_s[0]
  p1_a2 = p1Marker_s[1]
  p2_a1 = p2Marker_s[0]
  p2_a2 = p2Marker_s[1]
  
  
  
  if p1_a1 != p1_a2:
    if p2_a1 != p2_a2:
      if (p1_a1 == p2_a1 and p1_a2 == p2_a2) or (p1_a2 == p2_a1 and p1_a1 == p2_a2): #ab x ab
	selectedType = 'B3.7'
	if (p1_a1 == p2_a1 and p1_a2 == p2_a2):
	  a = p1_a1
	  b = p1_a2
	else:
	  a = p1_a2
	  b = p1_a1
	
      elif p1_a1 != p2_a1 and p1_a2 != p2_a2 and p1_a2 != p2_a1 and p1_a1 != p2_a2: #ab x cd
	selectedType = 'A.1'
	a = p1_a1
	b = p1_a2
	c = p2_a1
	d = p2_a2
	
      elif (p1_a1 == p2_a1 and p1_a2 != p2_a2) or (p1_a1 == p2_a2 and p1_a2 != p2_a1): #ab x ac
	selectedType = 'A.2'
	if (p1_a1 == p2_a1 and p1_a2 != p2_a2):
	  a = p1_a1
	  b = p1_a2
	  c = p2_a2
	else:
	  a = p1_a1
	  b = p1_a2
	  c = p2_a1
	
      else:
	selectedType = 'F' #selection failed
	print('DeterminMarkerType failed within (p1_a1 != p1_a2) and (p2_a1 != p2_a2): ' + p1_a1 + p1_a2 + p2_a1 + p2_a2)
	
    else: #p2_a1 == p2_a2
      if p2_a1 == '.': #ab x o
	selectedType = 'D1.11'
	a = p1_a1
	b = p1_a2	
	
      elif p1_a1 != p2_a1 and p1_a2 != p2_a1: #ab x c
	selectedType = DifferentiateBetweenA3_D19(allMarkerData,p1_a1,p1_a2,p2_a1)
	a = p1_a1
	b = p1_a2
	c = p2_a1
	
      elif p1_a1 == p2_a1 or p1_a2 == p2_a1:#ab x a
	if p1_a1 == p2_a1:
	  selectedType = DifferentiateBetweenB15_D110(allMarkerData,p1_a1,p1_a2)    
	  a = p1_a1
	  b = p1_a2
	else:
	  selectedType = DifferentiateBetweenB15_D110(allMarkerData,p1_a2,p1_a1)    
	  a = p1_a2
	  b = p1_a1
	  
      else:
	selectedType = 'F' #selection failed
	print('DeterminMarkerType failed within (p1_a1 != p1_a2) and (p2_a1 == p2_a2): ' + p1_a1 + p1_a2 + p2_a1 + p2_a2)
	
  else: #p1_a1 == p1_a2
    if p2_a1 != p2_a2:
      if p1_a1 == '.': #o x ab
	selectedType = 'D2.16'
	a = p2_a1
	b = p2_a2
	
      elif p1_a1 == p2_a1 or p1_a1 == p2_a2: # a x ab
	if p1_a1 == p2_a1:
	  selectedType = DifferentiateBetweenB26_D215(allMarkerData,p1_a1,p2_a2)
	  a = p2_a1
	  b = p2_a2
	  
	else:
	  selectedType = DifferentiateBetweenB26_D215(allMarkerData,p1_a1,p2_a1)
	  a = p2_a2
	  b = p2_a1
	  
      elif p1_a1 != p2_a1 and p1_a1 != p2_a2: #c x ab
	selectedType = 'D2.14'
	a = p2_a1
	b = p2_a2
	c = p1_a1	
      
      else:
	selectedType = 'F' #selection failed
	print('DeterminMarkerType failed within (p1_a1 == p1_a2) and (p2_a1 != p2_a2): ' + p1_a1 + p1_a2 + p2_a1 + p2_a2)       
	
    else: #p2_a1 == p2_a2
      if p2_a1 == '.': #a x o
	selectedType = 'D1.13'
	a = p1_a1
	
      elif p1_a1 == '.': #o x a
	selectedType = 'D2.18'
	a = p2_a1
	
      elif p1_a1 != p2_a1: #a x b or b x a
	selectedType = DifferentiateBetweenA4_D112_D217(allMarkerData,p1_a1,p2_a1)
	if selectedType == 'D1.12':
	  a = p2_a1
	  b = p1_a1
	else: 
	  a = p1_a1
	  b = p2_a1
	
      elif p1_a1 == p2_a1: #a x a 
	selectedType = 'C.8'      
	a = p1_a1
	
      else:
	selectedType = 'F' #selection failed
	print('DeterminMarkerType failed within (p1_a1 == p1_a2) and (p2_a1 == p2_a2): ' + p1_a1 + p1_a2 + p2_a1 + p2_a2)    

  return {'Type' : selectedType, 'a': a, 'b':b , 'c':c, 'd':d}
  
#Differentiates between A.3 and D1.9 marker types based on all marker data
#ab x c cross
def DifferentiateBetweenA3_D19(markerData,a,b,c):

  ac_test = False
  a_test = False
  bc_test = False
  b_test = False
  

  for marker in markerData:
    marker_s = marker.split('/')
    if marker_s[0] == a and marker_s[1] == c:
      ac_test = True
    if marker_s[0] == a and marker_s[1] == a:
      a_test = True
    if marker_s[0] == b and marker_s[1] == c:
      bc_test = True
    if marker_s[0] == b and marker_s[1] == b:
      b_test = True
      
  if ac_test and a_test and bc_test and b_test:
    return 'A.3'
  elif ac_test and b_test and c_test and not bc_test:
    return 'D1.9'
  else:
    print('DifferentiateBetweenA3_D19 failed')   
    print(markerData)
    return 'F'
   
#Differentiates between B1.5 and D1.10 marker types based on all marker data
#ab x a cross
def DifferentiateBetweenB15_D110(markerData,a,b):
  
  ab_test = False
  a_test = False
  b_test = False
  
  for marker in markerData:
    marker_s = marker.split('/')
    if marker_s[0] == a and marker_s[1] == b:
      ab_test = True
      
    if marker_s[0] == a and marker_s[1] == a:
      a_test = True

    if marker_s[0] == b and marker_s[1] == b:
      b_test = True
    
    if ab_test and a_test and b_test:
      return 'B1.5'
    elif ab_test and a_test and not b_test:
      return 'D1.10'
    else:
      print('DifferentiateBetweenB15_D110 failed')   
      print(markerData)
      return 'F'
    
#Differentiates between B2.6 and D2.15 marker types based on all marker data
#a x ab cross
def DifferentiateBetweenB26_D215(markerData,a,b):  
    
  ab_test = False
  a_test = False
  b_test = False
  
  for marker in markerData:
    marker_s = marker.split('/')
    if marker_s[0] == a and marker_s[1] == b:
      ab_test = True
      
    if marker_s[0] == a and marker_s[1] == a:
      a_test = True

    if marker_s[0] == b and marker_s[1] == b:
      b_test = True
    
    if ab_test and a_test and b_test:
      return 'B2.6'
    elif ab_test and a_test and not b_test:
      return 'D2.15'
    else:
      print('DifferentiateBetweenB26_D215 failed')   
      print(markerData)
      return 'F'
  
#Differentiates between A.4, D1.12 and D2.17 marker types based on all marker data
#axb or bxa cross
def DifferentiateBetweenA4_D112_D217(markerData,a,b):    
  
  ab_test = False
  a_test = False
  b_test = False
  o_test = False
  
  
  for marker in markerData:
    marker_s = marker.split('/')
    if marker_s[0] == a and marker_s[1] == b:
      ab_test = True
      
    if marker_s[0] == a and marker_s[1] == a:
      a_test = True

    if marker_s[0] == b and marker_s[1] == b:
      b_test = True
      
    if marker_s[0] == '.' or marker_s[1] == '.':
      o_test = True
      
  if ab_test and a_test and b_test and o_test:
    return 'A.4'
  #Next two tests check for axb or bxa by checking a_test and b_test
  elif ab_test and not a_test and b_test and not o_test:
    return 'D1.12'
  elif ab_test and a_test and not b_test and not o_test:
    return 'D2.17'
  else:
    print('DifferentiateBetweenA4_D112_D217 failed')   
    print(markerData)
    return 'F'
  
#marker determination is done based on patterens within the sibling markers.
def DeterminMarkerType_unknownParent(markerDict):
  
  #returned information
  selectedType = ''
  a = ''
  b = ''
  c = ''
  d = ''
  
  
  #band type counters
  c_a = 0
  c_b = 0
  c_c = 0
  c_d = 0
  c_ab = 0
  c_ac = 0
  c_ad = 0
  c_bc = 0
  c_bd = 0
  

  #run though all marker data and fill out the marker types
  for marker in markerDict:
    split_marker = marker.split('/')
    
    #if missing marker continue
    if split_marker[0] == '.':
      continue
    
    #if the data is already set then continue
    if split_marker[0] in [a,b,c,d] and split_marker[1] in [a,b,c,d]:
      continue
    
    
    if split_marker[0] == split_marker[1]:
      if a == '':
	a = split_marker[0]
      elif b == '':
	b = split_marker[0]
      elif c == '':
	c = split_marker[0]
      elif d == '':
	d = split_marker[0]
      else:
	print('DeterminMarkerType_unknownParent: Too many band types found. Failed.')
	return {'Type' : 'F', 'a': a, 'b':b , 'c':c, 'd':d}
   
    else:
      
      if split_marker[0] in [a,b,c,d]:
	if a == '':
	  a = split_marker[1]
	elif b == '':
	  b = split_marker[1]
	elif c == '':
	  c = split_marker[1]
	elif d == '':
	  d = split_marker[1]
	else:
	  print('DeterminMarkerType_unknownParent: Too many band types found. Failed.')
	  return {'Type' : 'F', 'a': a, 'b':b , 'c':c, 'd':d}
      elif split_marker[1] in [a,b,c,d]:
	if a == '':
	  a = split_marker[0]
	elif b == '':
	  b = split_marker[0]
	elif c == '':
	  c = split_marker[0]
	elif d == '':
	  d = split_marker[0]
	else:
	  print('DeterminMarkerType_unknownParent: Too many band types found. Failed.')
	  return {'Type' : 'F', 'a': a, 'b':b , 'c':c, 'd':d}
      else:
	if a == '' and b == '' and c == '' and d == '':
	  a = split_marker[0]
	  b = split_marker[1]
	elif a != '' and b == '' and c == '' and d =='':
	  b = split_marker[0]
	  c = split_marker[1]
	elif a != '' and b != '' and c == '' and d =='':
	  c = split_marker[0]
	  d = split_marker[1]
	else:
	  print('DeterminMarkerType_unknownParent: Too many band types found. Failed.')
	  return {'Type' : 'F', 'a': a, 'b':b , 'c':c, 'd':d}
	
  #count each time a specific band appears in the marker data
  for marker in markerDict:
    
    split_marker = marker.split('/')
    
    #if missing marker continue
    if split_marker[0] == '.':
      continue
    
    
    if split_marker[0] == split_marker[1]:
      if split_marker[0] == a:
	c_a += 1
      elif split_marker[0] == b:
	c_b += 1
      elif split_marker[0] == c:
	c_c += 1
      elif split_marker[0] == d:
	c_d += 1
      else:
	print('DeterminMarkerType_unknownParent: Unknown marker encountered during count')
	return {'Type' : 'F', 'a': a, 'b':b , 'c':c, 'd':d}
    else:
      if (split_marker[0] == a and split_marker[1] == b) or (split_marker[1] == a and split_marker[0] == b):
	c_ab += 1
      elif (split_marker[0] == a and split_marker[1] == c) or (split_marker[1] == a and split_marker[0] == c):
	c_ac += 1
      elif (split_marker[0] == a and split_marker[1] == d) or (split_marker[1] == a and split_marker[0] == d):
	c_ad += 1
      elif (split_marker[0] == b and split_marker[1] == c) or (split_marker[1] == b and split_marker[0] == c):
	c_bc += 1
      elif (split_marker[0] == b and split_marker[1] == d) or (split_marker[1] == b and split_marker[0] == d):
	c_bd += 1
      else:
	print('DeterminMarkerType_unknownParent: Unknown marker combo encountered during count')
	print(marker + ' ' + a + ' ' + b + ' ' + c + ' ' + d )
	return {'Type' : 'F', 'a': a, 'b':b , 'c':c, 'd':d}
    
  
  #use count and type specification to determin marker type
  if a != '' and b == '': #type C or D1.13 or D2.18
    selectedType = 'C'
  elif a != '' and b != '' and c == '': #type A.4 or B1.5 or B2.6 or B3.7 or D1.10 or D1.11 or D1.12 or D2.15 or D2.16 or D2.17
    
    if c_a > 0 and c_b > 0 and c_ab == 0: #D1.11 or D2.16
      selectedType = 'D1.11'
    elif (c_a > 0 and c_b == 0 and c_ab > 0): #D1.10 or D1.12 or D1.15 or D1.17
      selectedType = 'D1.10'
    elif (c_a == 0 and c_b > 0 and c_ab > 0): #D1.10 or D1.12 or D1.15 or D1.17 as well
      selectedType = 'D1.10'
      t_a = a
      a = b
      b = t_a
    elif c_a/2 >= c_b and c_a/2 >= c_ab: #B1.5 or B2.6
      selectedType = 'B1.5'
    elif c_ab/2 >= c_a and c_ab/2 >= c_b: #B3.7
      selectedType = 'B3.7'
    elif c_a > 0 and c_b > 0 and c_ab > 0: #A.4
      selectedType = 'A.4'
    else:	
      print('DeterminMarkerType_unknownParent: a + b selection failed.')
      print(str(c_a) + ' ' + str(c_b) + ' ' + str(c_c) + ' ' + str(c_ab) + ' ' + str(c_ac) + ' ' + str(c_ad) + ' ' + str(c_bc) + ' ' + str(c_bd))
      selectedType = 'F'
  elif a != '' and b != '' and c != '' and d == '': #type A.2 or A.3 or D1.9 or D2.14
    if c_a == 0 and c_b == 0 and c_c == 0 and c_ab == 0 and c_ac > 0 and c_bc > 0 : #D1.9 or D2.14
      selectedType = 'D1.9'
    elif c_a > 0 and c_b == 0 and c_c == 0 and c_ab > 0 and c_ac > 0 and c_bc > 0 : #A.2
      selectedType = 'A.2'
    elif c_a > 0 and c_b > 0 and c_c == 0 and c_ab == 0 and c_ac > 0 and c_bc > 0 : #A.3 
      selectedType = 'A.3'
    else:
      print('DeterminMarkerType_unknownParent: a + b + c selection failed.')
      selectedType = 'F'
  elif a != '' and b != '' and c != '' and d != '': #type A.1
    selectedType = 'A.1'
  
  
  return {'Type' : selectedType, 'a': a, 'b':b , 'c':c, 'd':d}

  
  
def GenerateMarkerEntry(markerName,markerTypeDict,markerData,p1_index,p2_index):
  
  entryString = ''
  
  entryString += '*' + markerName + ' ' + markerTypeDict['Type'] + '   '
  
  
  for mkrIdx in range(0,len(markerData)):
    #skip parent markers
    if mkrIdx == p1_index or mkrIdx == p2_index:
      continue
    
    marker_s = markerData[mkrIdx].split('/')
    
    if marker_s[0] == markerTypeDict['a'] and marker_s[1] == markerTypeDict['c']:
      entryString += 'ac,'
    elif marker_s[0] == markerTypeDict['a'] and marker_s[1] == markerTypeDict['d']:
      entryString += 'ad,'
    elif marker_s[0] == markerTypeDict['b'] and marker_s[1] == markerTypeDict['c']:
      entryString += 'bc,'
    elif marker_s[0] == markerTypeDict['b'] and marker_s[1] == markerTypeDict['d']:
      entryString += 'bd,'
    elif marker_s[0] == markerTypeDict['a'] and marker_s[1] == markerTypeDict['a']:
      entryString += 'a,'
    elif marker_s[0] == markerTypeDict['b'] and marker_s[1] == markerTypeDict['a']:
      entryString += 'ab,'
    elif marker_s[0] == markerTypeDict['a'] and marker_s[1] == markerTypeDict['b']:
      entryString += 'ab,'
    elif marker_s[0] == markerTypeDict['b'] and marker_s[1] == markerTypeDict['b']:
      entryString += 'b,'
    elif marker_s[0] == '.' and marker_s[1] == '.':
      entryString += '-,'
    else:
      entryString += '-,'
      print('GenerateMarkerEntry encountered unexpected marker: ' + markerData[mkrIdx] + ' for ' + markerName + ' (' + markerTypeDict['a'] + markerTypeDict['b'] + markerTypeDict['c'] + markerTypeDict['d'] + ')' )
      
  #remove last , in string
  entryString = entryString[:-1]

  #verify entryString is valid
  validBandList = list()
  
  if markerTypeDict['Type'] == 'A.1':
    validBandList.extend(list(['ac','ad','bc','bd']))
  elif markerTypeDict['Type'] == 'A.2':
    validBandList.extend(list(['a','ac','ba','bc']))
  elif markerTypeDict['Type'] == 'A.3':
    validBandList.extend(list(['ac','a','bc','b']))
  elif markerTypeDict['Type'] == 'A.4':
    validBandList.extend(list(['ab','a','b','o']))
  elif markerTypeDict['Type'] == 'B1.5' or markerTypeDict['Type'] == 'B2.6':
    validBandList.extend(list(['ab','a','b']))
  elif markerTypeDict['Type'] == 'B3.7':
    validBandList.extend(list(['a','ab','b']))
  elif markerTypeDict['Type'] == 'C.8':
    validBandList.extend(list(['a']))
  elif markerTypeDict['Type'] == 'D1.9':
    validBandList.extend(list(['ac','bc']))
  elif markerTypeDict['Type'] == 'D1.10':
    validBandList.extend(list(['a','ab']))
  elif markerTypeDict['Type'] == 'D1.11':
    validBandList.extend(list(['a','b']))
  elif markerTypeDict['Type'] == 'D1.12':
    validBandList.extend(list(['ab','a']))
  elif markerTypeDict['Type'] == 'D1.13':
    validBandList.extend(list(['a','o']))
  elif markerTypeDict['Type'] == 'D2.14':
    validBandList.extend(list(['ac','bc']))
  elif markerTypeDict['Type'] == 'D2.15':
    validBandList.extend(list(['a','ab']))
  elif markerTypeDict['Type'] == 'D2.16':
    validBandList.extend(list(['a','b']))
  elif markerTypeDict['Type'] == 'D2.17':
    validBandList.extend(list(['ab','a']))
  elif markerTypeDict['Type'] == 'D2.18':
    validBandList.extend(list(['a','o']))
  
  StartIndex =  entryString.index(',')
  
  if entryString[StartIndex - 2] == ' ':
    StartIndex -= 1
  else:
    StartIndex -= 2
  
  entryStringBandSplit = entryString[StartIndex:].split(',')
  
  for testPoint in entryStringBandSplit:
    if testPoint in validBandList or testPoint == '-':
      continue
    else:
      print('GenerateMarkerEntry - Unexpected marker in entryString: ' + testPoint + ' ' + markerTypeDict['Type'] + ' ' + markerName)
      print(entryString)
      return ''
      
      
  return entryString
    
     
    
  
  
  



  
  
  
  