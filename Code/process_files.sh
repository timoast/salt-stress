#! /bin/bash

# move each file to it's own folder
# move each folder to scratch
# remove duplicates
# index
# move back to working_data

alias picard='/home/tstuart/working_data/Tools/jre1.8.0_77/bin/java -jar /home/tstuart/working_data/Tools/picard-tools-2.1.1/picard.jar'

wd=`pwd`

for myfile in ./*.bam; do
    fname=(${myfile//.bam/ })
    mkdir $fname
    mv $myfile $fname
    mv $fname /home/tstuart/scratch
    cd /home/tstuart/scratch/$fname
    picard MarkDuplicates I=$myfile O=$fname.rmdup.bam M=dup_metrics.txt \
	   REMOVE_DUPLICATES=true
    samtools index $fname.rmdup.bam
    cd ..
    mv $fname $wd
    cd $wd
done

