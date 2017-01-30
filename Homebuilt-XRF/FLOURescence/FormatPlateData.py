##Written by Daniel Sanche
##Modified to change output directory on 10/27/2016
##Modified to clean up file output code 1/23/2017

import numpy as np
import pandas as pd
from os import listdir, walk, mkdir
from os.path import isfile, join, isdir, basename
import sys

#extracts the Ydata (Ydata, Yfit, energy) from a single .fit data file
def extractYdata(file):
    energyArr = []
    yfitArr = []
    ydatArr = []
    reading = False
    for line in file:
        row_delim = " = "
        if not reading and line.rstrip() == "[result]":
            # start reading under the [result] section
            reading = True
        elif reading and len(line.split(row_delim)) == 1:
            #stop reading when you reach the next section
            break
        elif reading:
            #look for yfit, energy, and ydata
            title, data = line.rstrip().split(row_delim);
            #strip brackets
            data = data.replace("]", "").replace("[", "").split()
            if title == "energy":
                energyArr = data
            elif title == "ydata":
                ydatArr = data
            elif title == "yfit":
                yfitArr = data
    combinedArr = np.column_stack((energyArr, ydatArr, yfitArr))
    df = pd.DataFrame(combinedArr, columns=["Energy", "Ydata", "Yfit"])
    return df

#extracts the peak fit from a single .fit data file
def extractPeakFit(file, spectraNum = None, key_txt="K.KL3]"):
    row_delim = " = "
    reading = False
    elementArr = []
    dataArr = []
    for line in file:
        if "]" in line and "[" in line:
            #we found a header. Check if it's one we're interested in
            if key_txt in line:
                reading = True
                title = line.rstrip()
                element = title.split()[0].replace("[result.", "")
                elementArr += [element]
            else:
                reading = False
        elif reading and len(line.split(row_delim)) == 2:
            # look for data
            title, data = line.rstrip().split(row_delim);
            if title == "fitarea":
                # strip brackets
                data = data.replace("]", "").replace("[", "")
                dataArr += [data]
    combinedArr = np.column_stack((elementArr, dataArr))
    spectraLabel = "Spectra"
    if spectraNum is not None:
        spectraLabel = spectraLabel + " " + str(spectraNum)
    df =pd.DataFrame(combinedArr, columns=["Element Label", spectraLabel])
    return df

#takes the yData dataframes from multiple .fit files and combines them into a single dataframe
def combineYdata(ydatDict):
    spectra = ydatDict.keys()
    spectra = sorted(spectra)
    key_titles = []
    dataframes = []
    for i in spectra:
        spectraTitle = "Spectra " + str(i)
        key_titles += [spectraTitle]
        df = ydatDict[i]
        dataframes += [df]
    combinedDf = pd.concat(dataframes, axis=1, keys=key_titles)
    return combinedDf

#takes the peakfit dataframes from multiple .fit files and combines them into a single dataframe
def combinePeakData(peakdatDict):
    spectra = peakdatDict.keys()
    spectra = sorted(spectra)
    first = True
    for i in spectra:
        df = peakdatDict[i]
        if first:
            first = False
            combinedData = df
        else:
            combinedData = pd.merge(combinedData, df, on="Element Label")
    return combinedData

#iterates though all .fit files in a directory, merges their data into a single dataframe, and saves it to disk
def ExtractDirectory(file_dir, output_dir, delim=","):
    dirsDict = {}
    ydatdir = join(output_dir,'spectra')
    peakfitdir = join(output_dir, 'peakfits')

    if not isdir(ydatdir):
        mkdir(ydatdir)
    if not isdir(peakfitdir):
        mkdir(peakfitdir)

    all_files = [f for f in listdir(file_dir) if isfile(join(file_dir, f))]
    for root, dirs, files in walk(file_dir):
        for file in files:
            if ".fit" in file:
                if root in dirsDict:
                    ydatDict = dirsDict[root][0]
                    peakfitDict = dirsDict[root][1]
                else:
                    ydatDict = {}
                    peakfitDict = {}
                spectra_num = int(file.replace(".fit", "").split(".")[-1])
                file = open(join(root, file))
                yDat = extractYdata(file)
                ydatDict[spectra_num] = yDat
                peakFit = extractPeakFit(file, spectra_num)
                peakfitDict[spectra_num] = peakFit
                dirsDict[root] = [ydatDict, peakfitDict]
    
    for root in dirsDict.keys():
        ydatDict = dirsDict[root][0]
        peakfitDict = dirsDict[root][1]
        if len(ydatDict) > 0 and len(peakfitDict) > 0:
            ydatCombined = combineYdata(ydatDict)
            ydat_filename = join(ydatadir, (basename(root).replace('_trans.csv_FITDIR', '') + "_ydata.csv"))
            ydatCombined.to_csv(ydat_filename, sep=delim, index=False)
            print('Writing...' + ydat_filename)
            peakcombined = combinePeakData(peakfitDict)
            peakfit_filename = join(ydatadir, (basename(root).replace('_trans.csv_FITDIR', '') + "_peakfit.csv"))
            peakcombined.to_csv(peakfit_filename, sep=delim, index=False)
            print('Writing...' + peakfit_filename)