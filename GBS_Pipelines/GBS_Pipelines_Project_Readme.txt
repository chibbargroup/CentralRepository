GBS Pipelines Project


This project is a collection of scripts that perform various GBS (Genotype by Sequencing) tasks on a set of 
illumina style DNA reads. The Documents folder contains various papers and manuals that were used in developing the
scripts.  


BT_ST_Pipeline.sh 

  This is a bash script that uses bowtie2 and samtools to align the illumina reads to a specified refrence genome
  then uses samtools to format and convert the results to a SAM alignment format.  This script is only designed
  for unpaired reads.  Paired reads will need a different script.
  
  Usage:
    There are three setting within the script that must be adjusted before the analysis can be run correctly:
    
      BOWTIE_GENOME_LOC="./velvet_assembly/HV_TestRef_cutoffAuto" #Bowtie refrence created by bowtie2-build
      SAMTOOLS_GENOME_LOC="./velvet_assembly/contigs.fa" #Samtools refrence created by sametools faidx
      NUM_THEADS="8" #Number of threads to allow analysis to use.  Used by bowtie and trimmomatic
    
    The script it self it run by a normal command line call with the name of the fastq read file as a argurment.
    In order to automate the process the 'find' command can be used in the following manner:
      
	find -iname '*.fastq' -exec ./BT_ST_Pipeline.sh {} \;
	
    The files that are passed into the script need to be already trimmed of the barcodes and seperated into different files based
    on the sample.
    
UNEAK_Pipeline.sh

  This is a bash script that uses the Tassle UNEAK pipline to perform refrenceless GBS analysis.  It requires 4 input arguments
  to run:
  
    ./UNEAK_Pipeline <Input Files> <Work Directory> <Barcode File> <Enzyime>
    
  The work Directory must be setup prior to running the pipeline script.  More information can be found in the UNEAK manuals
  
OneMap_File_Generator_f2.py

  This file is used to generate the input files for the R linkage analysis package OneMap.  It takes a tab version of a VCF file as an
  input and outputs a file that can be fed into the Onemap_Linkage.R script. The script needs to be run within python by calling the main()
  function.  This version creates the files for f2 crosses, (backcross, RIL, f2...)

  
OneMap_File_Generator_Outcross.py

  This file is used to generate the input files for the R linkage analysis package OneMap.  It takes a tab version of a VCF file as an
  input and outputs a file that can be fed into the Onemap_Linkage.R script. The script needs to be run within python by calling the main()
  function.  This version creates the files for outcrosses.
  
  
Onemap_Linkage.R

  This is a R script that performs linkage analysis on GBS data. It uses the OneMap package and follows the examples in OneMap manual. 
  The script takes unix style inputs with the following options:
  
    --lod <value>, Sets the LOD value for the analysis (Default: 3)
    --maxrf <value>, Sets the Max recombonational frequency (Default: 0.5)
    --test, The script will stop before the full mapping is done.  Used to check the number and size of groups.
    --f2 <String>, Input file in f2 format, cannot be used with --oc
    --oc <String>, Input File in outcross format, cannot be used with --f2
    
    
RNA-Seq_Pipeline_PE.sh

    This is a script for running RNA-Seq analysis on paired end data. The setup of the files is similar to normal GBS but we use TopHat instead of bowtie to do the alignemnts and cufflinks to do the transcriptome finding.  This script requires an extra refrence list file for use.  The file needs to be in the following format:
    
        <Name of Refrence> <location of refrence files minus the .fa> <location of the gff files for the refrence>
    
        A /mnt/data1/refrence_genomes/Triticum_Aestivum/Triticum_aestivum.IWGSC1.0+popseq.30.dna.chromosome.A /mnt/data1/refrence_genomes/Triticum_Aestivum/Triticum_aestivum.IWGSC1.0+popseq.30.chromosome.A.gff3
        B /mnt/data1/refrence_genomes/Triticum_Aestivum/Triticum_aestivum.IWGSC1.0+popseq.30.dna.chromosome.B /mnt/data1/refrence_genomes/Triticum_Aestivum/Triticum_aestivum.IWGSC1.0+popseq.30.chromosome.B.gff3
        D /mnt/data1/refrence_genomes/Triticum_Aestivum/Triticum_aestivum.IWGSC1.0+popseq.30.dna.chromosome.D /mnt/data1/refrence_genomes/Triticum_Aestivum/Triticum_aestivum.IWGSC1.0+popseq.30.chromosome.D.gff3
        
    GENOME_LOC_FILE variable needs to be setup before scripts is called and then its called with the paired-end files as inputs.
    
    See the RNA-Seq_AnalysisMethod in documents folder for more detail.
        
  
  
  