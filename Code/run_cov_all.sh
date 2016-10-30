#! /bin/bash

for myfile in ./*; do
    if [ -d $myfile ]; then
	cd $myfile
	python ../../Code/get_cov.py -f $myfile.rmdup.bam -o $myfile.tsv -s 100000
	cd ..
    fi
done

