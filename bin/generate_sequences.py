#!/usr/bin/env python

#################################################
# Generates sequences from different modules 	#
#  with potentially random linker sequences		#
#################################################

import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument('--xrRNA_file', 	'-x', type=str,	help='Input file with xrRNA motifs')
parser.add_argument('--spacer_file', 	'-s', type=str,	help='Input file with spacer sequences')
parser.add_argument('--IRES_file', 	'-i', type=str,	help='Input file with IRES sequences')
parser.add_argument('--CDS_file', 	'-c', type=str,	help='Input file with CDS sequences')
parser.add_argument('--output', 	'-o', type=str,	help='Output file')
parser.add_argument('--length','-l', type=int, default=30, help='Length of 5\' overhang ')
args = parser.parse_args()

########################
########################
###    parameters    ###
########################
########################

base_illegal_sequences = [
		"ATTTAAT",
		"GAAGAC",
		"GGTCTC",
		"CGTCTC",
		"CACNNNGTG",
		"TAATACGACTCACTATAG",
		"GCTAGTTATTGCTCAGCGG",
		"ATG",
		"TGA",
		"TAG",
		"TAA",
		"AAGMAMAACA",
		"CAGAAACACT",
		"CTCAAGCAACT",
		"AAACNAAKA",
		"AGTANATTNTSAKSACT",
		"CTCWMCWAGWAMTTT",
		"CAAWAYAWACGRATCTCAA",
		"AATMACWAMTCWMAMSACAACA",
		"CAGCRACRCAWMAGSAWNAAC",
		"TGTCSTKATYAARGCNNAATCASNAACA"
	]

##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main(sequence_file,output,length,stem_file,linker_file):
	errx('hi')



#######################
#######################
###    functions    ###
#######################
#######################
# Print error message and exit script	
def errx(message):
    print(message)
    exit(1)

##########################
### starts main script ###
##########################
main(sequence_file=args.sequence_file,
	output=args.output,
	length=args.length,
	stem_file=args.stem_file,
	linker_file=args.linker_file)


