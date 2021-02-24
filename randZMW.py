#!/usr/bin/env python

#--------------------------------------------------------------------------------------------------------------
# Downsampling PacBio reads by ZMW
# Author: Cameron Watson
# Last Updated: 23 Feb 2021
#--------------------------------------------------------------------------------------------------------------

import argparse
import re
import numpy as np
import gzip
import subprocess

#--------------------------------------------------------------------------------------------------------------
# USER INPUT
#--------------------------------------------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(description = "A tool for downsampling Pacbio reads by random ZMW in FASTA format",
    add_help = True)

    parser.add_argument("-f", "--fasta", 
    help = "A FASTA file of PacBio reads with the ZMW denoted in each header line", required = True)

    return parser.parse_args()

args = get_args()

#--------------------------------------------------------------------------------------------------------------
# FUNCTIONS
#--------------------------------------------------------------------------------------------------------------

class fastaRecord:

    __slots__ = ['header', 'zmw', 'seq']
    
    def __init__(self, record):
        '''
        instantiate fastaRecord objects w/ fasta header line, ZMW number, sequence
        '''
        self.header = record[0]
        self.zmw = self.header.split("/")[1]
        self.seq = record[1]

    def outputSelected(self, zmw_list, output_fh):
        '''
        write out fasta record if zmw is in the randomly selected list
        '''
        if self.zmw in zmw_list:
            output_fh.write(self.header)
            output_fh.write(self.seq)

#--------------------------------------------------------------------------------------------------------------
# MAIN
#--------------------------------------------------------------------------------------------------------------

# bash subprocess, grep for fasta headers
headers = subprocess.Popen(['grep','^>', args.fasta], 
    shell = False, stderr = subprocess.PIPE, stdout = subprocess.PIPE, text = True)

# pipe grep stdout to cut, cut for ZMW number
zmws = subprocess.Popen(['cut','-d', '/', '-f', '2'], stdin = headers.stdout, 
    shell = False, stderr = subprocess.PIPE, stdout = subprocess.PIPE, text = True)

# taking a random 75% of ZMWs
zmws = np.array(zmws.communicate()[0].split()) # array of all ZMWs

zmws = np.unique(zmws) # only unique ZMW numbers

subSize = int(len(zmws) * 0.75) # randomly select 75% of ZMWs

final_zmw = np.random.choice(zmws, size = subSize, replace = False)

# open input and output files
fh = open(args.fasta, "r")

basename = re.split("/", args.fasta)[-1]

prefix = re.split('\.', basename)[0]

outName = "./" + prefix + "_downsampled.fa.gz"

outFile = gzip.open(outName, "wt")

# iterate through fasta file 
ln = 0
current_record = ["",""]
for line in fh:
    if ln == 0: # captures first line
        current_record[0] = line
    elif re.match("^>", line): # captures all subsequent header lines
        current_record = fastaRecord(current_record) # see fastaRecord class
        current_record.outputSelected(final_zmw, outFile)
        current_record = ["",""]
        current_record[0] = line
    else: # sequence lines
        current_record[1] += line
    ln += 1

# handle EOF
current_record = fastaRecord(current_record) 
current_record.outputSelected(final_zmw, outFile)


fh.close()
outFile.close()