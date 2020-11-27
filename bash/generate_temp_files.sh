#!/usr/bin/env bash

# assign command line arguments to variables used in this script; $1 is first argument, $2 is second, and $3 is the third argument specified
file=$1
dir=$2
chrom=$3
# check if previous temporary file is present; if so, delete
if [ -f $dir/chrom_* ]
then
    rm $dir/chrom_*
fi
# sort the reads in the current chromosome by strand
samtools view -F 0x10 $dir/$file $chrom > $dir/chrom_${chrom}.txt
samtools view -f 0x10 $dir/$file $chrom >> $dir/chrom_${chrom}.txt
