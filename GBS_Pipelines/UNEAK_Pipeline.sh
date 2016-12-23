#!/bin/bash
#GBS Pipeline Using UNEAK in Tassle 3
#http://www.maizegenetics.net/#!tassel/c17q9
#Written by Craig Irvine (cri646@mail.usask.ca)
#21 April 15


echo "UNEAK GBS Pipeline Version 0.1"
echo "Writen by Craig Irvine"
echo "21 April 15"
echo ""
echo "Start Time:"
date

#Check number of arguments
if [ $# -ne 4 ];
  then echo "Invalid Number of Arguments"
  exit 1
fi


#settings
TASSLE_RUN_PIPE_LOC="/usr/bin/tassel3-standalone/run_pipeline.pl"
FASTQC_LOC="/usr/bin/FastQC/fastqc"
INPUT_SEQ_FILES_LOC=$1
INPUT_WORK_DIR=$2
INPUT_KEY_FILE_LOC=$3
ENZIME=$4

echo "Creating Work Directories"
#Create working directory if it doesn't exist and copy over files.
if [ ! -d $INPUT_WORK_DIR ];
  then mkdir $INPUT_WORK_DIR
fi

#Create Output Log Directory
if [ ! -d $INPUT_WORK_DIR/Output_logs ];
  then mkdir $INPUT_WORK_DIR/Output_logs
fi


#Create all working directories
$TASSLE_RUN_PIPE_LOC -fork1 -UCreatWorkingDirPlugin -w $INPUT_WORK_DIR -endPlugin -runfork1 > $INPUT_WORK_DIR/Output_logs/UCreatWorkingDirPlugin.log

#copy the key file to the key directory
cp $INPUT_KEY_FILE_LOC $INPUT_WORK_DIR/key/


echo "Trimming Sequence Files"
#trim the input sequence files and place results in /Illumina/ directory
#all files need to be in *_fasta.gz format
for path in $INPUT_SEQ_FILES_LOC/*.gz; do
    echo $path
    java -jar /usr/share/java/trimmomatic-0.32.jar SE -threads 8 -phred33 $path $INPUT_WORK_DIR/Illumina/$(basename $path) ILLUMINACLIP:TruSeq3-SE.fa:2:30:6 SLIDINGWINDOW:4:15 MINLEN:36 CROP:170 2> $INPUT_WORK_DIR/Output_logs/trim_output_$(basename $path).log
done

echo "Trimming done"

echo "FastQC"
$FASTQC_LOC $INPUT_WORK_DIR/Illumina/*.gz 
echo "FastQC Done"



echo "Running UNEAK Pipline"
#run rest of UNEAK pipeline on working directory
$TASSLE_RUN_PIPE_LOC -fork1 -UFastqToTagCountPlugin -w $INPUT_WORK_DIR -e $ENZIME -endPlugin -runfork1 > $INPUT_WORK_DIR/Output_logs/UFastqToTagCountPlugin.log

$TASSLE_RUN_PIPE_LOC -fork1 -UMergeTaxaTagCountPlugin -w $INPUT_WORK_DIR -c 5 -x 60000000 -m 120000000 -endPlugin -runfork1 > $INPUT_WORK_DIR/Output_logs/UMergeTaxaTagCountPlugin.log

$TASSLE_RUN_PIPE_LOC -fork1 -UTagCountToTagPairPlugin -w $INPUT_WORK_DIR -e 0.03 -endPlugin -runfork1 > $INPUT_WORK_DIR/Output_logs/UTagCountToTagPairPlugin.log

$TASSLE_RUN_PIPE_LOC -fork1 -UTagPairToTBTPlugin -w $INPUT_WORK_DIR -endPlugin -runfork1 > $INPUT_WORK_DIR/Output_logs/UTagPairToTBTPlugin.log

$TASSLE_RUN_PIPE_LOC -fork1 -UTBTToMapInfoPlugin -w $INPUT_WORK_DIR -endPlugin -runfork1 > $INPUT_WORK_DIR/Output_logs/UTBTToMapInfoPlugin.log

$TASSLE_RUN_PIPE_LOC -fork1 -UMapInfoToHapMapPlugin -w $INPUT_WORK_DIR -mnMAF 0.05 -mxMAF 0.5 -mnC 0.8 -mxC 1 -endPlugin -runfork1 > $INPUT_WORK_DIR/Output_logs/UMapInfoToHapMapPlugin.log


echo "Analysis Done"
date
