#!/usr/bin/env python

#################################################
# Extracts sequences that somewhat match the	# 
#  wanted secondary structure					#
#################################################

import argparse
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument('--sequence_file', 	'-s', type=str,	help='Input sequence with predicted secondary structure')
parser.add_argument('--output', 	'-o', type=str,	help='Output file')
args = parser.parse_args()

########################
########################
###    parameters    ###
########################
########################

target_structure = '.....((((....))))..........((((...)))).....'
pattern_start = r'\.{1,10}'
pattern_first_hp = r'\({3,5}\.{4}\){3,5}'
pattern_stretch = r'\.{7,11}'
pattern_second_hp_stem_left = r'(?:\({3,7}|\({2,6}\.?\({1,5}|\({1,5}\.?\({2,6})'
pattern_second_hp_loop = r'\.{3,10}'
pattern_second_hp_stem_right = r'(?:\){3,7}|\){2,6}\.?\){1,5}|\){1,5}\.?\){2,6})'
pattern_second_hp = fr'{pattern_second_hp_stem_left}{pattern_second_hp_loop}{pattern_second_hp_stem_right}'
pattern_end = r'\.{1,10}'

pattern = fr'{pattern_start}{pattern_first_hp}{pattern_stretch}{pattern_second_hp}{pattern_end}'



##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main(sequence_file,output):
	# Read file
	with open(sequence_file, 'r') as file:
		lines = file.readlines()

	collected_entries = []

	# Read lines and generate a new entry for every 3 lines
	for i in range(0,len(lines),3):
		# read lines
		id = lines[i].rstrip()
		sequence_nuc = lines[i+1].rstrip()
		sequence_struct = lines[i+2].rstrip()

		# extract information 
		split_header 				= id.split()
		split_sequence_structure 	= sequence_struct.split()
		

		header_name 					= split_header[0][1:] if id.startswith('>') else errx('First line of each entry needs to start with >\nThis line is: ' + id)
		position_xrRNA	 				= split_header[1].split(':')[1]
		position_spacer				 	= split_header[2].split(':')[1]
		position_ires					= split_header[3].split(':')[1]
		position_cds					= split_header[4].split(':')[1]
		structure_current				= split_sequence_structure[0]
		# Remove parentheses from MFE
		mfe 							= split_sequence_structure[1].replace('(','').replace(')','')

		xrRNA_start		= int(position_xrRNA.split('-')[0])
		xrRNA_end		= int(position_xrRNA.split('-')[1])
		xrRNA_structure = structure_current[xrRNA_start:xrRNA_end + 1]

		matches = re.finditer(pattern,xrRNA_structure)

		for match in matches:
			print(f'Found match in {header_name}')
			print(f'Found match: {match.group()}')
			print(f'Startposition: {match.start()}')
			print(f'Endposition: {match.end()}')
			
			collected_entries = collected_entries.append(SeqRecord(Seq(sequence_nuc),
														id = f'{id} regex:{match.group()}'))

	if len(collected_entries) == 0:
		print('No matches found')
		exit(0)

	with open(output,'w') as output_handle:
		SeqIO.write(collected_entries, output_handle, "fasta")



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
	output=args.output)
