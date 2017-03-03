# FLOURescence Manual version 0.1

----

## Introduction and Scope

This manual is designed to help users who are new to the FLOUResence software package. The software is designed to be intuitive to use, but this guide will help users analyze XRF data through all steps of the analysis process.

## Some Initial Things - Software Limitations

The FLOUResence software has a few limitations which must be dealt with ##WHEN COLLECTING## the data. Mostly these limitations deal with file names. When naming files at the beamline, try to avoid the following:
* The plate name cannot contain the "\_" character
* Anything after the plate name (separated by an underscore) will not be read by the program; i.e., for plate1_something, the something will be ignored
* Replicates from a plate cannot have the same plate name; i.e., the program cannot differentiate between plate1_1 and plate1_2; instead rerun the plate using a name like plate1-2ndRep or something similar

**IT IS IMPORTANT TO PROPERLY LABEL YOUR PLATES AT THE BEAMLINE AS THIS WILL SAVE YOU A LOT OF WORK DURING THE ANALYSIS**

## Installing Required Software

The following software should be installed on your computer to run the analysis
* Transpose
* FLOUResence
* PyMCA

Windows installers for these files can be found in the Resources folder in the MolecCropQualCommon shared folder on the Jade drive.

When installing Transpose and FLOUResence, it is recommended that create desktop shortcuts to the files FLOUResence.exe and Transpose.exe found in their respective installation directories.

# The General Analysis Work Flow

There are several steps in the workflow for analyzing data from the VESPERS beamline (detailed descriptions provided below):
1. Collect and collate the data; put Io files and spectra files in separate folders
2. Create a sample labels csv file
3. Transpose the spectra data using the program Transpose
4. Fit the spectra data using PyMCA
5. Analyze the fits and generate reports using FLOUResence

Note that steps 4 and 5 have to be repeated twice; once to quantify Ca, the second time to quantify the rest of the elements

----

## Step 1: Put Io files and spectra files in separate folders

This part is relatively straightforward. Spectra files generally have the formate plate_SomeDetector_#.dat where SomeDetector describes the type of detector used to collect the data (e.g., FourElementVortex or 13ElGe) and # represents a number. Spectra files generally will __NOT__ include the phrases 'vline', 'raw', or 'hline' in the name.

The files containing the Io data (herein referred to as the Io files) are generally significantly smaller than the other files. Also, they will usually have the format Plate_#.dat.

## Step 2: Create a csv file containing the sample labels

This file must have the following format:
| Plate      | Spot 01     | Spot 02     | Spot 03     | ... | Spot 84     |
|:----------:|:-----------:|:-----------:|:-----------:|:---:|:-----------:|
| plate name | Sample Name | Sample Name | Sample Name | ... | Sample Name |


The row in grey must be included in the file. Where plate name is the name of the plate. (Note that the name here must __exactly match the plate name used in the spectra file__). The spots are in the order that the spectra were collected. Note that usually __spectra are collected by row, starting with the bottom row from the bottom__!

Empty spots must be labelled as : X, x, Empty, or empty
Calibration samples must be labelled as: Cal # (where # is the calibration point number, i.e., Cal 1, Cal 2, Cal 3, etc.); Note, this is __CASE SENSITIVE__, so __must be exactly Cal #__.

Finally, one last requirement: The file must be saved using the __COMMA SEPARATED VALUE (.csv) format__. In Excel, this option can be chosen from the list of formats.

Pro Tip: Try using a previous file containing sample labels as a template if you get stuck/are confused

## Step 3: Tranpose the spectra files

The data in the spectra files must be transposed so they can be read by PyMCA. To do this, run the program Transpose.

The program will ask for two inputs:
1. The folder that contains the spectra files
2. A location to save the transposed files

If the folder where you'd like to save the transposed files does not exist, the program will conviently make it for you.

Note that the save location cannot be the same as the folder containing the spectra files. 

Also note that you will have to put in the full path to the directory; i.e., for a folder placed on the desktop, the path would be:

> C:\Users\User_Name\Desktop\Folder

## Step 4: Fit the Spectra Using PyMCA

Please see PyMCA_Tutorial.doc for more information about this step. Briefly, you can use the Batch Fitting tool found under the tools menu to generate fits. Note that you will have to run this twice, once for the Ca analysis and once for the non-Ca analysis, as the parameters used for fitting each region are slightly different.

## Step 5: Analyze the fits using FLOUResence

Seeing as this is a manual for FLOUResence, this will be described as its own section.

----

# Running FLOUResence

## Starting the Program

The FLOUResence program can be run by simply double clicking on the created shortcut or double clicking the FLOUResence.exe file. This will open up a command prompt window which you will use to interact with the program.

At the start, the program will ask the user for several parameters. All inputs will require that the full path to the file or folder is input, i.e.:

> C:\Users\User_Name\Desktop\Folder

The inputs asked for are:
1. The directory containing the fit results generated by PyMCA (NOTE: This isn't the directory you chose when running step 4, but rather the folder called FIT that is contained within the previously chosen output directory)
2. The .csv file containing the sample label information
3. The directory containing the Io files
4. The directory where you would like to save the output results; if the output directory does not exist, the program will ask if you would like to create it

Once your parameter have been input, the program repeat them and ask if they are correct. Here you can type Y or N. If you type N, then you will be asked to re-input the parameters. Also, if you have already run the program, it will remember the input parameters and files used previously, and will ask if you'd like to use them. If you select Y, it will show the previous parameters, at which time you can change them by following the on-screen instructions.

Once the input paramters have been entered, the program will give you a number of options to choose from, as shown here:

> What would you like to do? Type in the number to indicate procedure: 
> 1) Format plate data
> 2) Normalize the data and add sample labels
> 3) Combine the normalized data
> 4) Perform the linear regresssion analysis
> 5) Calculate concentrations and perform error analysis
> 6) Generate plots of the fits and calibration curves
> 7) Do all of the above (this may overwrite existing data in the output directory)
> 8) Exit the program

Options can be selected by typing the desired option number and then pressing Enter.

The descirptions of each option will be enumerated below.

## Option 1: Format Plate Data

This step reads the data from fit files and collates the data into csv files for each plate. Two sets of files are created in this step: 1) A series of CSV files that contain the peak areas and 2) the XRF spectra and their fits, and these are saved to the directory. The files are saved in the following two directories:

> Output_Directory/peakfits

> Output_Directory/spectra

Please note that if any files are in these directories there is a chance they will be overwritten by the program.

## Option 2: Normalize the data and add sample labels

In this step the peak area files generated by step one are normalized and then the sample labels are added to the resulting file. Peak areas are normalized by dividing by the appropriate Io measurement (norm peak area = peakarea/io). (Io is the intensity of the incident X-ray beam.) Sample labels are then added, and the results are saved in the directory:

> Output_Directory/norm

## Option 3: Combine the normalized data

Up to this point, all the data files have been separated by plate. In this this step, these data are combined into 4 separate files:
1. calibration_data.csv - contains the data from the individual calibration standards
2. calibration_data_avg.csv - contains the averaged data from each calibration standard
3. sample_data.csv - contains all of the peak areas from the measured samples
4. sample_data_avg.csv - contains the averages from the replicates of each sample

These files are stored in the directory: 

> Output_Directory/combined_files

## Option 4: Perform the Linear Regression Analysis

This module generates the calibration curves for each element by performing a linear regression on the averaged calibration data (contained in the file calibration_data_avg.csv). 

During this step, the program will ask you for the concentration of each calibration sample. Simply input the concentration (no units) for each calibration point. If you have already input the concentrations, the program will remember and ask if you would like to use the previous values.

The calibration curve parameters and calibration concentrations will be saved in the directory:

> Output_Directory/calibration_results

Note that Limit of Detection (LOD) and Limit of Quantification (LOQ) values are contained in the folder 

## Option 5: Calculate concentrations and perform error analysis

The final step in the data processing loop, this uses calibration curve results to convert the the peak areas to analyte concentrations. The module also calculates the error associated with the measurement. Two files are saved during this process:

1. sample_concentrations.csv - contains the concentrations determined for each individual sample
2. avg_results.csv - contains the averaged results for each sample, also includes the error bar for each measurement

## Option 6: Generate plots of the fits and calibration curves

The plotting module generates plots of for XRF spectrum and fit. It also relabels the spectra files to include the sample names. The plots of the spectra are saved in the following directory:

> Output_Directory/spectra_plots

while the labelled spectra files are saved in the directory:

> Output_Directory/labeled_spectra

In addition to plotting the spectra, this module also generates plots of the calibration curves for each element studied. The calibration curve plots show the line of best fit and the averaged calibration data points. These plots are saved in the directory:

> Output_Directory/calibration_results/cal_curve_plots

__WARNING__: There are usually quite a few spectra (upwards of 1000), so this step can take a long time to run (usually it takes about 1 s/spectrum)

## Option 7: Do all of the above (this may overwrite existing data in the output directory)

This option will run steps 1-5 in one go. It currently does not run step 6.

## Option 8: Exit the program

Pretty self explanatory. Will close the program and the command window.

## Other Notes

Because this version of FLOUResence uses averaged calibration data rather than individual calibration data, it is possible that outliers can significantly affect the results. it is recommended that you check the combined calibration data to see if there are outliers. It may be necessary to manually change the info in the averaged calibration data file to ensure these outliers do not affect the results significantly. (Outliers for Fe seem to be relatively common, likely do to easy contamination issues.)
