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
#07 Oct 2015 - Changed command lines to support .fq.gz files instead of plain .fq file to save space
#10 Nov 2015 - Added ILLUMINACLIP option to trimmomatic to trim out the Illumina Atapters from the sequences
#             - Commented out vcf pileup of individual bam files to speed up processing time.
#31 Mar 2016 - Removed -trimlog option from trimmomatic to save space.
#            - Removed 170 bp crop from trimmomatic

echo "Bowtie/Samtools Pipeline Version 0.8"
echo "Writen by Craig Irvine"
echo ""


#settings
BOWTIE_GENOME_LOC="/mnt/data1/refrence_genomes/Triticum_Aestivum/TGACv1.30/Triticum_aestivum.TGACv1.30.dna.genome" #Bowtie refrence created by bowtie2-build
SAMTOOLS_GENOME_LOC="/mnt/data1/refrence_genomes/Triticum_Aestivum/TGACv1.30/Triticum_aestivum.TGACv1.30.dna.genome.fa" #Samtools refrence created by sametools faidx
NUM_THEADS="8" #Number of threads to allow analysis to use.  Used by bowtie and trimmomatic


#check that the input file exists
if [ ! -e $1 ]
then
  echo "No input file supplied"
  exit 1
fi

#check the refrences exist
#BOWTIE_GENOME_LOC
if [ ! -e "$BOWTIE_GENOME_LOC.1.bt2l" ]
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


echo "Processing input file $1";
date

#create output folder for input file
OUTPUT_FOLDER="$1_output"
mkdir $OUTPUT_FOLDER


#Trimming filters the reads by quality set the the SLIDINGWINDOW and MINLEN parts of the command.  
#Trimming requires fastQ file that contain the quality data, if fastA files are being used coment out the 
#trimmomatic command and uncomment the cp command.  This is to allow the chain of output file to not be broken
#and allow the rest of the script to continue to run.
echo "Trimming file $1"
java -jar /usr/share/java/trimmomatic-0.32.jar SE -threads $NUM_THEADS -phred33 $1 ./$OUTPUT_FOLDER/$1.trim.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:6 HEADCROP:20 SLIDINGWINDOW:4:15 MINLEN:36 2> ./$OUTPUT_FOLDER/trim_output.log
#cp $1 ./$OUTPUT_FOLDER/$1.trim
echo "Trimming Complete"


echo "Creating FastQC File of Trimmed Data"
/usr/bin/FastQC/fastqc -q ./$OUTPUT_FOLDER/$1.trim.gz
echo "FastQC done"


#Bowtie aligns the reads to the refrence sequence
# --end-to-end force a full read align, this can be changed to --local if partial aligns are alowed.
# If fastA files are used as input, a -f should be placed before the -U option.
echo "Bowtie $1"
bowtie2 --end-to-end -p $NUM_THEADS -x $BOWTIE_GENOME_LOC -U ./$OUTPUT_FOLDER/$1.trim.gz -S ./$OUTPUT_FOLDER/$1.trim.bowtie2.sam 2> ./$OUTPUT_FOLDER/bowtie2_output.log
echo "Bowtie complete"

#Samtools View converts the SAM output from bowtie2 to a BAM file.
echo "samtools view"
samtools view -bS ./$OUTPUT_FOLDER/$1.trim.bowtie2.sam > ./$OUTPUT_FOLDER/$1.trim.bowtie2.bam 2> ./$OUTPUT_FOLDER/samtools_view_output.log
echo "samtools view complete"

#Samtools sort... sorts the bam file by read name
echo "samtools sort"
samtools sort ./$OUTPUT_FOLDER/$1.trim.bowtie2.bam ./$OUTPUT_FOLDER/$1.trim.bowtie2.sort 2> ./$OUTPUT_FOLDER/samtools_sort_output.log
echo "samtools sort complete"

#Samtools rmdup removes duplicate read alignments.
echo "samtools rmdup"
samtools rmdup ./$OUTPUT_FOLDER/$1.trim.bowtie2.sort.bam ./$OUTPUT_FOLDER/$1.trim.bowtie2.sort.rmdup.bam 2> ./$OUTPUT_FOLDER/samtools_rmdup_output.log
echo "samtools rmdup complete"

#samtools index creates an index of the BAM file so it can be loaded into a viewer
echo "samtools index"
samtools index ./$OUTPUT_FOLDER/$1.trim.bowtie2.sort.rmdup.bam
echo "samtools index complete"

#Appends final BAM file to a list for latter use with samtools mpileup
echo "Appending bam file to list"
echo "./$OUTPUT_FOLDER/$1.trim.bowtie2.sort.rmdup.bam" >> bam_list.txt

#Creates a pileup for the single read file.
echo "samtools mpileup"
#samtools mpileup -Duf $SAMTOOLS_GENOME_LOC ./$OUTPUT_FOLDER/$1.trim.bowtie2.sort.rmdup.bam | bcftools view -vg - > ./$OUTPUT_FOLDER/$1.trim.bowtie2.sort.rmdup.raw.vcf
echo "samtools mpileup complete"

echo "Cleanup"
rm -f ./$OUTPUT_FOLDER/$1.trim.gz
rm -f ./$OUTPUT_FOLDER/$1.trim.bowtie2.sam
rm -f ./$OUTPUT_FOLDER/$1.trim.bowtie2.bam
rm -f ./$OUTPUT_FOLDER/$1.trim.bowtie2.sort.bam

echo "$1 Processing Complete"
date
echo ""
