#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

prefix=(${1//.fastq*/ })
cores=$3
index="/scratchfs/moliva/Bowtie2Index/genome"

# Trim the nextera adapters
sh /home/sbuckberry/working_data_01/bin/bbmap/bbduk2.sh \
in=$1 \
in2=$2 \
out=$prefix.tmp_R1.fq \
out2=$prefix.tmp_R2.fq \
rliteral=TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG,GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG ktrim=r \
mink=3 \
threads=$cores \
overwrite=true

## bowtie2 alignment
bowtie2 -q --no-unal --threads $cores -X 2000 \
    -x $index \
    -1 $prefix.tmp_R1.fq -2 $prefix.tmp_R2.fq \
    | samtools view -q 30 -bSu - \
    > $prefix.temp

samtools sort -@ $cores -O bam -T sorted $prefix.temp -o $prefix.tmp.bam

# Calculate coverage stats
# bash /home/sbuckberry/working_data_01/bin/bbmap/pileup.sh \
# in="$prefix".tmp.bam \
# out="$prefix".covstats \
# overwrite=t &

# Create the bam file index
samtools index $prefix.tmp.bam

#Filter the chrM data
samtools idxstats $prefix.tmp.bam | cut -f 1 | grep -vE 'Mt|Pt' | xargs samtools view -b $prefix.tmp.bam > $prefix.bam

# Get the insert size metrics and generate plots
samtools view -f66 $prefix.bam | cut -f 9 | sed 's/^-//' > $prefix.InsertSizeMetrics.txt

# call the plotting script
Rscript /home/sbuckberry/working_data_02/polo_project/polo_iPSC/ATACseq/analysis_scripts/plotInsertSizeHistogram.R $prefix.InsertSizeM\
etrics.txt $prefix

# Remove the tmp files
rm "$prefix".tmp_R1.fq "$prefix".tmp_R2.fq "$prefix".tmp.bam "$prefix".tmp.bam.bai \
"$prefix".InsertSizeMetrics.txt $prefix.temp
