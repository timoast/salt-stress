#!/bin/bash

# Added safety. If any command fails or a pipe breaks, the script will stop running.
set -eu -o pipefail -o verbose

bam="$1"
cores="$2"
prefix=$(echo "$bam" | sed 's/.bam//g')

# Sort bam file by read name
samtools sort -n -@ "$cores" -T temp -O bam "$prefix".rmdup.bam > "$prefix".sorted.bam

# Update/fix flags
samtools fixmate "$prefix".sorted.bam "$prefix".fixed.bam

# Create bedpe file
bedtools bamtobed -i "$prefix".fixed.bam -bedpe > "$prefix".tmp.bedpe

# Sort bedpe file, remove discordant reads, reads in pair more than 2k bp away from each other; add "chr" prefix
sortBed -i "$prefix".tmp.bedpe | awk -v OFS="\t" '{if ( !($1==$4) || $9==$10 || $5-$3>2000); else print "chr"$1,$2,$3,"chr"$4,$5,$6,$7,$8,$9,$10}' - > "$prefix".filter.bedpe

# Offset (+5 for + mapped reads, -4 for - mapped reads)
awk -v OFS="\t" '{if ($9=="+") print $1, $2+4, $3+4, $4, $5-5, $6-5, $7,$8,$9,$10; else print $1, $2-5,$3-5, $4, $5+4, $6+4, $7,$8,$9,$10}' "$prefix".filter.bedpe > "$prefix".bedpe

# Remove tmp files
rm "$prefix".sorted.bam "$prefix".fixed.bam "$prefix".tmp.bedpe "$prefix".filter.bedpe