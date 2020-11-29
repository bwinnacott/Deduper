#!/usr/bin/env python
import re
import os
import argparse
from collections import defaultdict
import subprocess
import random
import shutil

def possible_umis(umi_file):
    '''Function that takes in a filename consisting of a list of unique molecular indexes (UMIs) and generates a dictionary 
    from the list. This dictionary represents all the possible UMIs present in the experiment.'''
    # initiate dictionary to store UMIs; open UMI file
    umi_dict = {}
    umis = open(umi_file,'r')
    # iterate over lines in file and store UMIs in UMI_dict
    for line in umis:
        line = line.strip()
        umi_dict[line] = []
    # keys: UMIs; values: empty lists
    return umi_dict

def error_correct_umi(curr_umi,umi_dict):
    '''Function designed to perform error correction on a given UMI that is only 1 nt in Hamming distance away from 
    only 1 UMI contained within the list of known UMIs in UMI_dict. If more than one UMI in UMI_dict have a Hamming 
    distance of 1, then the UMI is not corrected.'''
    # initiate variable to store possible corrected UMI
    poss_umi = None
    # iterate over all known UMIs and initiate Hamming distance variable
    for umi in umi_dict:
        hamm_dist = 0
        # iterate over all nt in UMI of interest
        for i in range(len(curr_umi)):
            # check for mismatched nts; if more than two are present for current UMI in 
            # UMI_dict, then move to next UMI
            if curr_umi[i] != umi[i]:
                hamm_dist += 1
                if hamm_dist > 1:
                    break
        # if distance is 2
        if hamm_dist > 1:
            continue
        # if current UMI in UMI_dict has Hamming distance of 1, assign to variable
        if poss_umi is None:
            poss_umi = umi
        # if UMI is already assigned, no longer error correct
        else:
            raise KeyError
    # account for no change in variable or return the only UMI in UMI_dict with Hamming distance of 1
    if poss_umi is None:
        raise KeyError
    else:
        return poss_umi

def check_strand(flag):
    '''Function that takes in bitwise flag (integer) and returns whether current read mapped as reverse complement to reference.'''
    # set variable for reverse complementarity
    reverse_comp = False
    # check flag to see if current read was reverse complemented
    if ((int(flag) & 16) == 16):
        reverse_comp = True
    # return boolean
    return reverse_comp

def check_clipping(cigar,reverse=False):
    '''Function that takes in a read's cigar string and returns the number for which the reads leftmost mapping position should 
    be shifted, depending on the reference strand to which the read mapped.'''
    # enter if reverse=True
    if reverse:
        # store amount to shift leftmost mapping position; get cigar operations
        position_shift = 0
        cigar_reverse = re.findall(r'(\d+)(\w+?)',cigar)
        # iterate over cigar operations; continue if first operation is 'S' or if operation is 'I'; for all others ('M','D','N',last 'S')
        # add to variable
        for i,vals in enumerate(cigar_reverse):
            if i == 0 and vals[1] == 'S':
                continue
            elif vals[1] == 'I':
                continue
            else:
                position_shift += int(vals[0])
        # return value
        return position_shift
    # enter if reverse=False
    else:
        # get first 'S' operation if present and return value
        cigar_forward = re.search(r'[0-9]+S',cigar)
        if cigar_forward is None:
            return 0
        else:
            return int(cigar_forward.group()[:-1])

def check_position(left_pos,cigar,reverse=False):
    '''Function that takes in leftmost mapping position, cigar string, and boolean indicating the strand to which the read mapped and returns 
    the adjusted position.'''
    # if read is reverse strand, return leftmost mapping position added to value returned from 'check_clipping'
    if reverse:
        return int(left_pos) + check_clipping(cigar,reverse=True)
    # if forward strand, return leftmost mapping position subtracted by value returned from 'check_clipping'
    else:
        return int(left_pos) - check_clipping(cigar)

def convert_phred(letter):
    """Converts a single character (letter) into a phred score"""
    # converts the ASCII symbol to phred equivalent (-33)
    score = ord(letter) - 33
    
    return score

def choose_read(read_group,method):
    '''Function that takes in a list of reads and returns the correct read based on the method selected for choosing one out of many.'''
    # if the random method is chosen, enter if
    if method == 'random':
        # select a random index from the group and return the read corresponding to that index
        chosen_read_idx = random.randrange(0,len(read_group))
        chosen_read = read_group.pop(chosen_read_idx)
        
        return (chosen_read,read_group)
    # if quality method is chosen, enter if
    else:
        # initiate appropriate variables
        max_mean_score = 0
        idx = 0
        # iterate through all reads in group
        for ind,val in enumerate(read_group):
            # get total quality score by iterating through each character in quality string
            total = 0
            for char in val[10]:
                total += convert_phred(char)
            # get average score for current read
            avg = total/len(val[10])
            # if current average is greater than largest seen so far, swap out variables
            if avg > max_mean_score:
                max_mean_score = avg
                idx = ind
        # get read with highest average quality
        chosen_read = read_group.pop(idx)
        # return all reads; separate chosen read
        return (chosen_read,read_group)

def write_out_reads(read_groups,method,deduped_fh,duplicate_fh):
    '''Function to write out reads according to the method specified other than "first" for chosing a read from a group
    of duplicates.'''
    # gather stats to add to the overall count
    total_chosen = 0
    total_not_chosen = 0
    # iterate over groups of reads
    for _,reads in read_groups.items():
        # if group only has 1 read, write it out
        if len(reads) == 1:
            deduped_fh.write('\t'.join(reads[0]))
            total_chosen += 1
        # choose a read from a set of multiple; write out to appropriate file
        else:
            chosen,not_chosen = choose_read(reads,method)
            deduped_fh.write('\t'.join(chosen))
            total_chosen += 1
            for read in not_chosen:
                duplicate_fh.write('\t'.join(read))
                total_not_chosen += 1
    # return stats
    return total_chosen,total_not_chosen

def get_chroms(file):
    '''Function that takes in the input SAM file and returns the list of chromosomes for which the reads mapped to.'''
    # use subprocess package to run bash script 'get_chroms.sh'
    chroms = subprocess.run(['bash/get_chroms.sh',file,'deduper_temp'],stdout=subprocess.PIPE)
    chroms = chroms.stdout.decode('utf-8')
    chroms = chroms.split('\n')[:-1]
    # return the list of chromosomes
    return chroms

def add_sam_header(file,output_file):
    '''Function to take the header lines (starting with '@') from the input SAM file and place it at the beginning 
    of every output file generated from running this script.'''
    # open input file and iterate line by line
    with open(file,'r') as fh:
        for line in fh:
            # add current line to output file if it starts with '@'
            if line.startswith('@'):
                output_file.write(line)
            # exit the context manager
            else:
                return

def open_files(input_file):
    '''Function to open all files for writing reads out to. File handles are passed around.'''
    # initiate file handles for the three files generated from the program
    deduped_fh = open('{}_deduped.sam'.format(re.sub('.sam','',input_file)),'w')
    duplicate_fh = open('{}_duplicate.sam'.format(re.sub('.sam','',input_file)),'w')
    errors_fh = open('{}_UMIerrors.sam'.format(re.sub('.sam','',input_file)),'w')
    # for each output file, write out the SAM file header lines
    for output_fh in [deduped_fh,duplicate_fh,errors_fh]:
        add_sam_header(input_file,output_fh)
    # return file handles
    return deduped_fh,duplicate_fh,errors_fh

def deduplicate_reads(input_file,umi_file,deduped_fh,duplicate_fh,errors_fh,error_correct,method):
    '''Main function to deduplicate reads from an input SAM file. Also takes as arguments the UMI file specified at the 
    command line, the method for which to select a read from a group of duplicates, and the three file handles for 
    writing out to.'''
    # open the input SAM file and initiate some variables collecting program stats
    with open(input_file,'r') as fh:
        total_reads,num_deduped,num_duplicates,num_umierrors,num_error_corr = 0,0,0,0,0
        # create dictionary with list of known UMIs; variable denoting current strand; dictionary to hold read group info
        poss_umis = possible_umis(umi_file)
        reverse = False
        # second dictionary to hold read info; keys: (UMI,adj position) tuples, vals: list of lists containing read content
        read_groups = defaultdict(list)
        # iterate over lines in the file
        for line in fh:
            total_reads += 1
            current_read = line.split('\t')
            # check to see if now dealing with reverse reads mapped to current chromosome
            if check_strand(current_read[1]) and not reverse:
                reverse = True
                # reset dictionaries
                poss_umis = possible_umis(umi_file)
                # write out all remaining reads in the 'read_groups' dict
                if method:
                    total_dedup,total_duplic = write_out_reads(read_groups,method,deduped_fh,duplicate_fh)
                    num_deduped += total_dedup
                    num_duplicates += total_duplic
                # reset dictionary to capture reads mapping to reverse strand
                read_groups = defaultdict(list)
            # get current UMI and check if it is known
            current_umi = current_read[0].split(':')[-1]
            if current_umi not in poss_umis:
                # correct if specified at command line
                if error_correct:
                # if unknown, try to error correct
                    try:
                        current_umi = error_correct_umi(current_umi,poss_umis)
                        num_error_corr += 1
                    # if error correction fails, write out to file containing UMI error reads; increment stat and move to next read
                    except:
                        errors_fh.write('\t'.join(current_read))
                        num_umierrors += 1
                        continue
                else:
                    errors_fh.write('\t'.join(current_read))
                    num_umierrors += 1
                    continue
            # get the adjusted position of current read
            adjusted_position = check_position(current_read[3],current_read[5],reverse=reverse)
            # check if adjusted position is in dictionary
            if adjusted_position in poss_umis[current_umi]:
                # if user wants first read chosen from group of duplicates, write out as duplicate
                if not method:
                    duplicate_fh.write('\t'.join(current_read))
                    num_duplicates += 1
                # if other method selected, record read information
                else:
                    read_groups[(current_umi,adjusted_position)].append(current_read)
            # add the current adjusted position to UMI dict at appropriate key
            else:
                poss_umis[current_umi].append(adjusted_position)
                # if user wants first read chosen from group of duplicates, write out a deduped read
                if not method:
                    deduped_fh.write('\t'.join(current_read))
                    num_deduped += 1
                # if other method selected, record read information
                else:
                    read_groups[(current_umi,adjusted_position)].append(current_read)
        # this code chunk accounts for the set of reverse reads stored after the entire file is parsed through
        if method:
            total_dedup,total_duplic = write_out_reads(read_groups,method,deduped_fh,duplicate_fh)
            num_deduped += total_dedup
            num_duplicates += total_duplic
    # return the stats
    return total_reads,num_deduped,num_duplicates,num_umierrors,num_error_corr

def generate_stats(total_reads,num_deduped,num_duplicates,num_umierrors,num_error_corr):
    '''Function that takes in all stats recorded throughout this program and spits out a report.'''
    # open the output file and write stats out
    with open('deduper_stats.txt','w') as f:
        f.write('Total number of reads processes: ' + str(total_reads) + '\n')
        f.write('Number of deduplicated reads: ' + str(num_deduped) + '\n')
        f.write('Percent of total reads that are deduplicated: ' + str(round((num_deduped/total_reads)*100,3)) + '%' + '\n')
        f.write('Number of duplicate reads: ' + str(num_duplicates) + '\n')
        f.write('Percent of total reads that are duplicates: ' + str(round((num_duplicates/total_reads)*100,3)) + '%' + '\n')
        f.write('Number of reads with UMI errors: ' + str(num_umierrors) + '\n')
        f.write('Percent of total reads that have UMI errors: ' + str(round((num_umierrors/total_reads)*100,3)) + '%' + '\n')
        f.write('Number of UMI corrected reads: ' + str(num_error_corr) + '\n')
        f.write('Percent of total reads that are UMI corrected: ' + str(round((num_error_corr/total_reads)*100,3)) + '%' + '\n')

def main():
    '''Main call to run the program.'''
    # create a temporary directory to store intermediate files generated and used throughout the program
    if not os.path.isdir('./deduper_temp'):
        os.mkdir('./deduper_temp')
    # initiate command line arguments/options
    parser = argparse.ArgumentParser(description='Tool used to remove PCR duplicates from a SAM file. Three separate files are produced which \
                                    include: 1) a file containing deduplicated reads, 2) a file containing the duplicated reads, and 3) a file \
                                    containing reads for which the UMI was not identified due to sequencing errors.')
    parser.add_argument('-f','--file',action='store',required=True,type=str,help='Specifies the input SAM file for which removal of PCR duplicates is desired.')
    parser.add_argument('-p','--paired',action='store_true',help='Designates the input file consists of paired end reads.')
    parser.add_argument('-u','--umi',action='store',type=str,help='Designates file containing list of UMIs specific to the reads contained within SAM file.')
    parser.add_argument('-e','--error_correct_umis',action='store_true',help='Error correction of UMIs is performed when specified.')
    parser.add_argument('-m','--method',action='store',choices=['random','quality'],type=str,help='Specifies the method used to determine \
                        the read kept from a group of duplicates. Method "random" chooses a random read from a group of duplicates and method "quality" \
                        selects the read with the highest average quality of a group of reads. If not specified, first duplicate encountered is kept.')
    # get the arguments specified at the command line
    args = parser.parse_args()
    # account for program's inability to handle paired end data or randomers (i.e., runs without known UMIs)
    if args.paired:
        parser.error('Current build of the program cannot handle paired end data processing. Rerun with single end data.')
    if not args.umi:
        parser.error('Current build of the program requires a file containing a list of known UMIs. Rerun with file of known UMIs specified with -u flag.')
    # run the program; first open files for writing out to
    deduped_fh,duplicate_fh,errors_fh = open_files(args.file)
    total_reads,num_deduped,num_duplicates,num_umierrors,num_error_corr = 0,0,0,0,0
    # iterate over the list of chromosomes and create temporary files sorted by strand
    print('...generating bam files...')
    for x in get_chroms(args.file):
        print('...starting chromosome {}...'.format(x))
        subprocess.run(['bash/generate_temp_files.sh','bam_input_sorted.bam','deduper_temp',x],stdout=subprocess.PIPE)
        # deduplicate the reads in each temp file and gather stats
        total,deduped,duplicates,umierrors,corrected = deduplicate_reads('deduper_temp/chrom_{}.txt'.format(x),args.umi,deduped_fh,duplicate_fh, \
            errors_fh,args.error_correct_umis,args.method)
        total_reads += total
        num_deduped += deduped
        num_duplicates += duplicates
        num_umierrors += umierrors
        num_error_corr += corrected
        print('...finished deduplicating chromosome {}.'.format(x))
    # create stats report
    generate_stats(total_reads,num_deduped,num_duplicates,num_umierrors,num_error_corr)
    # close all output files
    for f in [deduped_fh,duplicate_fh,errors_fh]:
        f.close()
    # remove the temporary directory created in the beginning of the program
    shutil.rmtree('./deduper_temp')

if __name__ == "__main__":
    main()