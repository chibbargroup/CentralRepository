#!/bin/bash
#Written by Craig Irvine (cri646@mail.usask.ca)
#23 Feb 2016

#Change log
# 23 Feb 2016: - Inital Version for RNA-Seq pipline.

echo "RNA-Seq Pipeline for PE Version 0.1"
echo "Writen by Craig Irvine"
echo ""


#settings
#Genomes must have the same name for tophat
GENOME_LOC_FILE="./RefList"
NUM_THEADS="8" #Number of threads to allow analysis to use.  Used by bowtie and trimmomatic


#check that the input files exists
if [ ! -e $1 ]
then
  echo "No input files supplied"
  exit 1
fi

if [ ! -e $2 ]
then
  echo "Paired Input File Missing"
  exit 1
fi


echo "Processing paired input files $1 $2";
date

#determining common string between input files
COMMON_INPUT="$(printf "%s\\n%s\\n" "$1" "$2" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/')"
echo "Common Input String: $COMMON_INPUT" 

#create output folder for input file
OUTPUT_FOLDER="$COMMON_INPUT"
mkdir $OUTPUT_FOLDER

#Trimming filters the reads by quality set the the SLIDINGWINDOW and MINLEN parts of the command.  
echo "Trimming files $1 $2"
java -jar /usr/share/java/trimmomatic-0.32.jar PE -threads $NUM_THEADS -phred33 $1 $2 ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1P.gz ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1U.gz ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2P.gz ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2U.gz HEADCROP:13 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36 2> ./$OUTPUT_FOLDER/trim_output.log
echo "Trimming Complete"

#FastQC
echo "Starting FASTQC"
/usr/bin/FastQC/fastqc -t $NUM_THEADS -o $OUTPUT_FOLDER $OUTPUT_FOLDER/*.gz
echo "Finished fastQC"

while read ref_line;
do
    echo "$ref_line"


    IFS=' '
    read -r -a arr <<< "$ref_line"

    OUTPUT_DIR=${arr[0]}
    GENOME_LOC=${arr[1]} 
    GFF_LOC=${arr[2]} 

    echo "Processing Refrences $GENOME_LOC"
    
    #create output folder for tophat
    OUTPUT_FOLDER="$COMMON_INPUT"
    mkdir $OUTPUT_FOLDER/$OUTPUT_DIR


    #check the refrences exist
    #BOWTIE_GENOME_LOC
    if [ ! -e "$GENOME_LOC.1.bt2" ]
    then
    echo "BOWTIE_GENOME_LOC refrence file missing"
    exit 1
    fi

    #SAMTOOLS_GENOME_LOC
    if [ ! -e "$GENOME_LOC.fa" ]
    then
    echo "SAMTOOLS_GENOME_LOC refrence file missing"
    exit 1
    fi
    
    #gff3 file
    if [ ! -e "$GFF_LOC" ]
    then
    echo "GFF3 File Missing"
    exit 1
    fi
    
    
    #Tophat - Run on paired sequences
    echo "Tophat $1 $2 Paired"
    tophat -p $NUM_THEADS -o $COMMON_INPUT/$OUTPUT_DIR -G $GFF_LOC $GENOME_LOC ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1P.gz ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2P.gz &> ./$OUTPUT_FOLDER/$OUTPUT_DIR/Tophat_output.log
    echo "Tophat Paired complete"
        
    
    #Cufflinks Start
    echo "Cufflink Assembly"
    cufflinks -q -p $NUM_THREDS -g $GFF_LOC -o $COMMON_INPUT/$OUTPUT_DIR  $COMMON_INPUT/$OUTPUT_DIR/accepted_hits.bam &> $COMMON_INPUT/$OUTPUT_DIR/cufflink_log
    echo "Cufflink Assembly Done"

    #Complie list of merged transcripts for each input genome
    echo "$COMMON_INPUT/$OUTPUT_DIR/merged.gtf" >> transcript_list_$OUTPUT_DIR


done<$GENOME_LOC_FILE







#echo "Cleanup"
#rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1P.gz
#rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1U.gz
#rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2P.gz
#rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2U.gz
#rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sam
#rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.bam
#rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.bam

echo "$1 $2 Processing Complete"
date
echo ""
