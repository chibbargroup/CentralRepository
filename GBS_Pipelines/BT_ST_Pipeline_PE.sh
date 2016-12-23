#!/bin/bash
#Written by Craig Irvine (cri646@mail.usask.ca)
#13 June 2014

#Change log
#16 June 2014 - Added time stamp at the begining and end of the script
#	      - Removed mpileup step and added creation of a list of all created .bam file
#		mpileup step is performed with all the .bam files at once.
#12 Dec 2014  - Added setting variables for Bowtie and samtools refrence locations
#	      - Added checks to ensure refrence files exist.
#	      - Commented out Trimming section for working on data that was already trimmed.
#21 Apr 2015 - Added setting variables for Number of Used Threads
#	     - Renamed script from NGS_Analysis.sh to BT_ST_Pipeline.sh	
#09 June 2015 - Documentation update
#	      - Added some error checking
#24 July 2015 - Rewrite to allow for paired end reads, will be treated as a seperate script from the SE version
#	      - Reset version number to 0.1 for PE script.
#23 Feb 2016  - Added FastQC analysis to pipeline.


echo "Bowtie/Samtools Pipeline for PE Version 0.2"
echo "Writen by Craig Irvine"
echo ""


#settings
BOWTIE_GENOME_LOC="/mnt/data1/refrence_genomes/Triticum_Aestivum/TA_D" #Bowtie refrence created by bowtie2-build
SAMTOOLS_GENOME_LOC="/mnt/data1/refrence_genomes/Triticum_Aestivum/Triticum_aestivum.IWGSC1.0+popseq.30.dna.chromosome.A.fa" #Samtools refrence created by sametools faidx
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

#check the refrences exist
#BOWTIE_GENOME_LOC
if [ ! -e "$BOWTIE_GENOME_LOC.1.bt2" ]
then
  echo "BOWTIE_GENOME_LOC refrence file missing"
  exit 1
fi

#SAMTOOLS_GENOME_LOC
if [ ! -e "$SAMTOOLS_GENOME_LOC" ]
then
  echo "SAMTOOLS_GENOME_LOC refrence file missing"
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
java -jar /usr/share/java/trimmomatic-0.32.jar PE -threads $NUM_THEADS -phred33 $1 $2 ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1P.gz ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1U.gz ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2P.gz ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2U.gz HEADCROP:10 SLIDINGWINDOW:4:20 MINLEN:36 2> ./$OUTPUT_FOLDER/trim_output.log
echo "Trimming Complete"


echo "FastQC Analsysis"
/usr/bin/FastQC/fastqc -q -t $NUM_THEADS ./$OUTPUT_FOLDER/*.gz
echo "FastQC Finished"


#Bowtie aligns the reads to the refrence sequence
# --end-to-end force a full read align, this can be changed to --local if partial aligns are alowed.
# If fastA files are used as input, a -f should be placed before the -U option.
echo "Bowtie $1 $2"
bowtie2 --end-to-end -p $NUM_THEADS -x $BOWTIE_GENOME_LOC -1 ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1P.gz -2  ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2P.gz -S ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sam 2> ./$OUTPUT_FOLDER/bowtie2_output.log
echo "Bowtie complete"

#Samtools View converts the SAM output from bowtie2 to a BAM file.
echo "samtools view"
samtools view -bS ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sam > ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.bam 2> ./$OUTPUT_FOLDER/samtools_view_output.log
echo "samtools view complete"

#Samtools sort... sorts the bam file by read name
echo "samtools sort"
samtools sort ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.bam ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort 2> ./$OUTPUT_FOLDER/samtools_sort_output.log
echo "samtools sort complete"

#Samtools rmdup removes duplicate read alignments.
echo "samtools rmdup"
samtools rmdup ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.bam ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.rmdup.bam 2> ./$OUTPUT_FOLDER/samtools_rmdup_output.log
echo "samtools rmdup complete"

#samtools index creates an index of the BAM file so it can be loaded into a viewer
echo "samtools index"
samtools index ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.rmdup.bam
echo "samtools index complete"

#Appends final BAM file to a list for latter use with samtools mpileup
echo "Appending bam file to list"
echo "./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.rmdup.bam" >> bam_list.txt

#Creates a pileup for the single read file.
#echo "samtools mpileup"
#samtools mpileup -Duf $SAMTOOLS_GENOME_LOC ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.rmdup.bam | bcftools view -vg - > ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.rmdup.raw.vcf
#echo "samtools mpileup complete"

echo "Cleanup"
rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1P.gz
rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_1U.gz
rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2P.gz
rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim_2U.gz
rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sam
rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.bam
rm -f ./$OUTPUT_FOLDER/$COMMON_INPUT.trim.bowtie2.sort.bam

echo "$1 $2 Processing Complete"
date
echo ""
