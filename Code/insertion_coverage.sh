#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

bedpe="$1"
prefix=$(echo "$bedpe" | sed 's/.bedpe//g')

# Create bed file from bedpe

cut -f 1-2,6-9 "$bedpe" > "$prefix".bed


# Create 300pb reads centered on each Tn5 insertion event (5' and 3' of paired-end reads)

awk -v OFS="\t" '{if ($2>50) print $1,$2-150,$2+149 RS $1,$3-149,$3+150; else print $1,$3-149,$3+150}' "$prefix".bed > "$prefix"_insertion.tmp.bed

bedClip "$prefix"_insertion.tmp.bed /home/moliva/working_data/genomes/TAIR10/chromosomesize.txt "$prefix"_insertion.fixed.tmp.bed

sort -k1,1 -k2,2n -k3,3n "$prefix"_insertion.fixed.tmp.bed > "$prefix"_insertion.bed


# Calculate coverage for Tn5 insertions

scalingFactor=$(cat "$prefix"_insertion.bed | awk 'END {print 1000000/NR}')

bedtools genomecov -bg \
-scale "$scalingFactor" \
-i "$prefix"_insertion.bed \
-g /home/moliva/working_data/genomes/TAIR10/chromosomesize.txt > "$prefix"_insertion.bedgraph

# Create BigWig file

bedGraphToBigWig "$prefix"_insertion.bedgraph /home/moliva/working_data/genomes/TAIR10/chromosomesize.txt "$prefix"_insertion.bigwig




rm "$prefix"_insertion.tmp.bed "$prefix"_insertion.fixed.tmp.bed

