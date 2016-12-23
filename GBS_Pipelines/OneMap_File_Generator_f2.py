#OneMap_File_Generator_f2.py
#This script will generate from a collection of vcf SNP calls a input OneMap_File_Generator
#for OneMap in order to facilitate linkage mapping. For F2, Backcrosses, and RIL Analysis.
#Input: A file containing a list of file locations of the files to be analyized.
#Output: Onemap Input file in the proper format.
#Version 0.2

#10Nov15: Ver 0.2 - Added Change log
#                 - Added new function for generating RIL marker entries
#                 - Ensured newlines were stripped from file names before appending to fileList
#                 - Added new funcion to count missing data and remove entries if a threshold is violated.



#main method for the script

def main(inputVCFFile,crossType,inputParent1Loc = '',inputParent2Loc = ''):
  import os.path
  
  #check that file exits
  if os.path.isfile(inputVCFFile) == False:
    print('Specified File does not exist: ' + inputVCFFile)
    return
  
  #check that the crossType is valid
  if not(crossType == 'f2 backcross' or crossType == 'f2 intercross' or crossType == 'ri self' or crossType == 'ri sib'):
    print('Specified Cross Type is invalid, must be f2 backcross, f2 intercross, ri self or ri sib')
    return
  
  
  
  
  #initialize data var
  markerDict = {}
  statsDict = {}
  fileList = list()   
   
  LoadVCFTabFile(inputVCFFile,markerDict,statsDict,fileList)    
    
  #if inputParent1Loc adn inputParent2Loc are not set, use first two samples in fileList
  if inputParent1Loc == '':
    inputParent1Loc = fileList[0]
    
  if inputParent2Loc == '':
    inputParent2Loc = fileList[1]
    
    
  p1_index = fileList.index(inputParent1Loc)
  p2_index = fileList.index(inputParent2Loc)  
 
  print('P1 Index ' + str(p1_index) + ' P2 Index ' + str(p2_index))
    

  numOfMarkers = 0  
  outFileString = ''
  
  for key in markerDict:
    if crossType == 'f2 backcross':
      markerEntry = GenerateBackcrossMarkerEntry(key,markerDict[key])
      numOfMarkers += 1
      numOfSamples = len(fileList) 
    elif crossType == 'f2 intercross':
      markerEntry = GenerateIntercrossMarkerEntry(key,markerDict[key],p1_index,p2_index)
      numOfMarkers += 1
      numOfSamples = len(fileList)-2
    elif crossType == 'ri self':
      if markerDict[key][p1_index] == markerDict[key][p2_index]:
          print('Parent SNPs match, skiping marker ' + markerDict[key][p1_index] + ' ' + markerDict[key][p2_index] + ' ' + key)
          continue
      
      p1Marker_s = markerDict[key][p1_index].split('/')
      p2Marker_s = markerDict[key][p2_index].split('/')
      
      if (p1Marker_s[0] != p1Marker_s[1]) or (p2Marker_s[0] != p2Marker_s[1]):
          print('Parent Markers are not homozygous, skipping ' + markerDict[key][p1_index] + ' ' + markerDict[key][p2_index] + ' ' + key)
          continue
      
      if (p1Marker_s[0] == '.') or (p2Marker_s[0] == '.'):
          print('Parent Missing Marker, skipping ' + markerDict[key][p1_index] + ' ' + markerDict[key][p2_index] + ' ' + key)
          continue
      
        
      markerEntry = GenerateRILMarkerEntry(key,markerDict[key],p1_index,p2_index)
      numOfMarkers += 1
      numOfSamples = len(fileList)
      
      if CountMissing(markerEntry,0.1,numOfSamples) == False:
          print('To many missing entries, removing marker:\n ' + markerEntry)
          numOfMarkers -= 1
          continue
      
    elif crossType == 'ri sib':
      markerEntry = GenerateBackcrossMarkerEntry(key,markerDict[key])
      numOfMarkers += 1
      numOfSamples = len(fileList)

    outFileString += markerEntry + '\n'
    
    


  
  #write out file
  outFile = open(inputVCFFile + '.onemap','w')
  
  outFile.write('data type ' + crossType + '\n')
  outFile.write(str(numOfSamples) + ' ' + str(numOfMarkers) + ' 0\n\n')
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
	  fileList.append(section.strip('\n'))
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
	  
  
   
#http://www.animalgenome.org/bioinfo/resources/manuals/carthagene/node15.html
    
  
def GenerateBackcrossMarkerEntry(markerName,markerData):
  
  entryString = ''
  
  entryString += '*' + markerName + ' '
  
  for marker in markerData:
    marker_s = marker.split('/')
    
    
    if marker_s[0] == '.':
      entryString += '-'   
    
    elif marker_s[0] == marker_s[1]:
      entryString += 'A'
      
    else:
      entryString += 'H'
      
  return entryString
    
  
def GenerateIntercrossMarkerEntry(markerName,markerData,p1_index,p2_index):
  
  entryString = ''
  
  entryString += '*' + markerName + ' '

  
  p1Marker_s = markerData[p1_idx].split('/')
  p2Marker_s = markerData[p2_idx].split('/')
  
  p1_allele = p1Marker_s[0]
  p2_allele = p2Marker_s[0]
  
  for markerIdx in range(0,len(markerData)):
    
    
    if markerIdx == p1_index or markerIdx == p2_index:
      continue
    
    marker_s = markerData[markerIdx].split('/')
    
    
    if marker_s[0] == '.':
      entryString += '-'   
    
    elif marker_s[0] == marker_s[1] and marker_s[0] == p1_allele:
      entryString += 'A'
      
    elif marker_s[0] == marker_s[1] and marker_s[0] == p2_allele:
      entryString += 'B'
      
    elif marker_s[0] == p1_allele and marker_s[1] == p2_allele:
      entryString += 'H'
      
    elif marker_s[0] == p1_allele:
      entryString += 'C'
      
    elif marker_s[0] == p2_allele:
      entryString += 'D'
      
    else:
      entryString += '-'
      
  return entryString
  
  
  
def GenerateRILMarkerEntry(markerName,markerData,p1_idx,p2_idx):
    
    
  entryString = ''
  
  entryString += '*' + markerName + ' '

  
  p1Marker_s = markerData[p1_idx].split('/')
  p2Marker_s = markerData[p2_idx].split('/')
  
  p1_allele = p1Marker_s[0]
  p2_allele = p2Marker_s[0]
  
  for markerIdx in range(0,len(markerData)):
    
    
    #if markerIdx == p1_idx or markerIdx == p2_idx:
     # continue
    
    marker_s = markerData[markerIdx].split('/')
    
    
    if marker_s[0] == '.':
      entryString += '- '   
    
    elif marker_s[0] == marker_s[1] and marker_s[0] == p1_allele:
      entryString += 'A '
      
    elif marker_s[0] == marker_s[1] and marker_s[0] == p2_allele:
      entryString += 'B '
      
    elif marker_s[0] == p1_allele and marker_s[1] == p2_allele:
      entryString += '- '
      
    elif marker_s[0] == p2_allele and marker_s[1] == p1_allele:
      entryString += '- '
      
    elif marker_s[0] == p1_allele and marker_s[1] != p1_allele:
      entryString += '- '
      
    elif marker_s[0] == p2_allele and marker_s[1] != p2_allele:
      entryString += '- '
      
    else:
      entryString += '- '
      
  return entryString
  
  
def CountMissing(entryString,max_missing,total_markers):
    
    missingCount = 0;
    
    for idx in range(0,len(entryString)):
        if entryString[idx] == '-':
            missingCount = missingCount + 1
            
            
    if missingCount > max_missing*total_markers:
        print('Missing Count: ' + str(missingCount) + '>' + str(max_missing*total_markers))
        return False
    
    else:
        #print('Missing Count: ' + str(missingCount) + '<' + str(max_missing*total_markers))
        return True
    
    
        