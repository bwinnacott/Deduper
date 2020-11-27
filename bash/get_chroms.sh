#!/usr/bin/env bash

# assign command line arguments to variables used in this script; $1 is first argument and $2 is second
file=$1
dir=$2
# create sorted bam files for efficient indexing
samtools view -bS $file > $dir/bam_input_file.bam
samtools sort $dir/bam_input_file.bam -o $dir/bam_input_sorted.bam
samtools index $dir/bam_input_sorted.bam
# get a list of all unique chromosomes associated with the input SAM file
cat $file | grep -v "^@" | cut -f 3 | sort -V | uniq
