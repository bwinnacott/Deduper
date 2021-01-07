# Deduper

## Overview
The goal of this tool is to perform reference-based PCR duplicate removal. For a given input 
SAM file, all PCR duplicates are removed so that only a single copy of each read is retained 
(default is first read encountered). Currently, the tool only accepts single-end, uniquely 
mapped reads as input. Additionally, a separate file containing a list of known UMIs is required 
to confirm the identity of duplicate reads (see *STL96.txt* for an example format). 

The program utilizes ```samtools``` to produce an intermediate SAM file sorted by mapping position, 
which is then broken up into subsequent temporary files by chromosome. Sets of reads mapped to 
a given chromosome are then sorted by strand to reduce the overall number of reads stored in 
memory. At the end of each chromosome, temporary files are deleted and created again in a temporary 
directory. 

## Requirements
- Python v3
- Samtools

## Output Files
Three output files are generated from this tool:
- SAM file containing **deduplicated reads**
- SAM file containing identified **PCR duplicates**
- SAM file containing **reads with sequencing errors** in UMI sequence

## Arguments
The following arguments are required for running the program:
- ```-f```, ```--file```: specifies absolute path to input SAM file for which removal of PCR 
duplicates is desired
- ```-u```, ```--umi```: designates absolute path to text file containing the list of known UMIs 
specific to the reads contained within the input SAM file

The following arguments are optional:
- ```-p```, ```--paired```: designates the input SAM file consists of paired end reads (**Note:** 
future option)
- ```-e```, ```--error_correct_umis```: error correction of UMI sequences is performed when 
specified
- ```-m```, ```--method```: Specifies the method used to determine the read kept from a group of 
duplicates. "random" chooses a read at random from a given group and "quality" selects the read 
with the highest average sequence quality (based on q-scores) from a group of reads.

