COUNTER=0
for d in */ ; do
    cd $d
    pwd
    qsub ../pileup_split_job.pbs
    cd ..
    COUNTER=$[$COUNTER +1]
done
