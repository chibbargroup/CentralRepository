#!/bin/bash

#This script starts multiple instances of the Structure command line application.
#It will run the application according to the settings on the files pointed to
#by the file locations.  
#max_processes should be set to the number of processors available for use.  
#see http://pritchardlab.stanford.edu/structure.html for STRUCTURE info

#Created 5-Aug-15 by Craig Irvine (cri646@mail.usask.ca)
#Version 0.1
#Change log
#
#

#File locations
structure_exe_loc='/usr/local/bin/Structure/console/structure'
mainparams_loc='./mainparams' #Needs to be setup for the specific input file
extraparams_loc='./extraparams' #No change needed usually
input_file='../Chickpea_setup/project_data' #Project data setup in Structure format
output_loc='./Results' #location where to place the results
log_loc='./Logs' #Location where to place the run logs.

#Settings
K_min=1
K_max=10
Iterations=10
max_processes=8
sleep_time='10m' #see Unix sleep command for possible settings. Current is 10 min


echo "Starting Structure Runs"
date


#check that log and results locations exists
if [ ! -e $output_loc ]
then
  echo "Creating Results Directory"
  mkdir -p $output_loc
fi

if [ ! -e $log_loc ]
then
  echo "Creating Log Directory"
  mkdir -p $log_loc
fi

command_list=() #Command list array to hold all generated structure runs


#Loop over all options and create the command line for the structure runs
for k in $(seq $K_min $K_max);
do
    for i in $(seq 1 $Iterations);
    do
    
    command_list+=("nohup $structure_exe_loc -K $k -m $mainparams_loc -e $extraparams_loc -i $input_file -o $output_loc/K$k-$i > $log_loc/K$k-$i.log &")
    
    done

done

#command and control loop, run though the command_list and run each command concurrently upto a max of max_processes
#monitor for when a process completes and run the next process in the list.
for i in "${command_list[@]}"
do

    #check if max_processes has been reached, if so wait 10 min and check again.
    
    num_processes="$(ps -ax | grep "$structure_exe_loc" | grep -v grep | wc -l)"
    while [ $num_processes -ge $max_processes ]
    do
        echo "Max processes reached ($num_processes), sleeping($sleep_time)"
        sleep $sleep_time
        num_processes="$(ps -ax | grep "$structure_exe_loc" | grep -v grep | wc -l)"
    done   
     
     
    echo "Running $i"
    eval $i
    sleep 1s  #Sleep for 1 seconds to allow processes to get upto speed before checking
    

done

date
echo "Structure Runs Complete"
