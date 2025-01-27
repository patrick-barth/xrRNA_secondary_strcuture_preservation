#!/usr/bin/env python

#################################################
# Splits an mRNA sequence to shorter fragments	#
#  depending on secondary structure. Longer 	# 
#  unpaired stretches are more likely to be 	#
#  used. Needs output file of RNAfold as input	#
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

length_linker = 4
length_stem = 6
length_loop = args.length_loop


##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main(sequence_file,output, length_loop):
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
		position_on_origin 				= split_header[1].split(':')[1]
		structure_of_origin_position 	= split_header[2].split(':')[1]
		stem_name 						= split_header[3].split(':')[1]
		linker_name						= split_header[4].split(':')[1]
		structure_current				= split_sequence_structure[0]
		# Remove parentheses from MFE
		mfe 							= split_sequence_structure[1].replace('(','').replace(')','')

		# Perform Checks if certain structure criteria are met
		
		# Check stem structure
		five_prime_stem_structure 	= structure_current[length_linker:(length_linker+length_stem)]
		three_prime_stem_structure	= structure_current[(len(structure_current)  - length_linker - length_stem):(len(structure_current) - length_linker)][::-1]
		stem_completely_paired		= False
		stem_paired_with_itself		= False
		stem_paired_nucleotides 	= None
		stem_unpaired_nucleotides 	= None

		# First Check if both stem sequences are fully paired
		if has_repeated_char(five_prime_stem_structure,'(',length_stem) and has_repeated_char(three_prime_stem_structure,')',length_stem):
			stem_completely_paired 		= True
			stem_paired_with_itself		= True
			stem_paired_nucleotides 	= length_stem
			stem_unpaired_nucleotides 	= 0
		# Check if stem is paired with itself
			# TODO: Check does not work perfectly. some instances are assigned as True even though the string is different
		elif not five_prime_stem_structure == three_prime_stem_structure[::-1].replace(')','('):
			stem_completely_paired		= False
			stem_paired_with_itself		= True
			stem_paired_nucleotides 	= five_prime_stem_structure.count('(')
			stem_unpaired_nucleotides 	= length_stem - stem_paired_nucleotides
		else:
			stem_completely_paired		= False
			stem_paired_with_itself		= False
			stem_paired_nucleotides 	= None
			stem_unpaired_nucleotides 	= None

		# Check linker structure
		five_prime_linker_structure = structure_current[0:length_linker]
		three_prime_linker_structure = structure_current[len(structure_current) - length_linker:] 
		linker_completely_unpaired = False
		linker_unpaired_nucleotides = None
		linker_paired_to_other_parts = False

		# Check if linker is completely free
		if has_repeated_char(five_prime_linker_structure,'.',length_linker) and has_repeated_char(three_prime_linker_structure,'.',length_linker):
			linker_completely_unpaired 		= True
			linker_unpaired_nucleotides 	= length_linker
			linker_paired_to_other_parts 	= False
		elif not five_prime_linker_structure == three_prime_linker_structure[::-1].replace(')','('):
			linker_completely_unpaired 		= False
			linker_unpaired_nucleotides 	= 'NA'
			linker_paired_to_other_parts 	= True
		else:
			linker_completely_unpaired 		= False
			linker_unpaired_nucleotides 	= five_prime_linker_structure.count('.')
			linker_paired_to_other_parts 	= False


		# Check loop structure
		loop_structure = structure_current[length_linker + length_stem:length_linker + length_stem + length_loop]
		loop_completely_unpaired 		= False
		loop_longest_unpaired_stretch 	= None
		loop_unpaired_nucleotides		= None
		loop_paired_nucleotides 		= None
		if has_repeated_char(loop_structure,'.', length_loop):
			loop_completely_unpaired 		= True
			loop_longest_unpaired_stretch 	= length_loop
			loop_unpaired_nucleotides		= length_loop
			loop_paired_nucleotides 		= 0
		else:
			loop_completely_unpaired 		= False
			loop_longest_unpaired_stretch 	= get_longest_consecutive_char_chain(structure_current,'.')
			loop_unpaired_nucleotides		= loop_structure.count('.')
			loop_paired_nucleotides 		= length_loop - loop_unpaired_nucleotides



		# Add entry
		collected_entries.append({
			'id':header_name,
			'position_on_origin':position_on_origin,
			'structure_of_origin_position':structure_of_origin_position,
			'stem_name':stem_name,
			'linker_name':linker_name,
			'structure_current':structure_current,
			'sequence_nuc':sequence_nuc,
			'mfe':mfe,
			'stem_completely_paired':stem_completely_paired,
			'stem_paired_nucleotides':stem_paired_nucleotides,
			'stem_unpaired_nucleotides':stem_unpaired_nucleotides,
			'stem_paired_with_itself':stem_paired_with_itself,
			'loop_completely_unpaired':loop_completely_unpaired,
			'loop_longest_unpaired_stretch':loop_longest_unpaired_stretch,
			'loop_unpaired_nucleotides':loop_unpaired_nucleotides,
			'loop_paired_nucleotides':loop_paired_nucleotides,
			'linker_completely_unpaired':linker_completely_unpaired,
			'linker_unpaired_nucleotides':linker_unpaired_nucleotides,
			'linker_paired_to_other_parts':linker_paired_to_other_parts
		})


	# write output
	# write header line
	output_string = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
		'id',
		'structure',
		'sequence',
		'stem_completely_paired',
		'stem_paired_with_itself',
		'stem_number_paired',
		'stem_number_unpaired',
		'loop_completely_unpaired',
		'loop_longest_unpaired_stretch',
		'loop_paired_nucleotides',
		'loop_unpaired_nucleotides',
		'linker_paired_to_other_parts',
		'linker_completely_unpaired',
		'linker_unpaired_nucleotides',
		'MFE',
		'loop_position_on_original_RNA',
		'local_structure_of_original_loop',
		'stem_name',
		'linker_name'
	)
	for entry in collected_entries:
		# add one line per entry
		output_string += '\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
			entry['id'],
			entry['structure_current'],
			entry['sequence_nuc'],
			entry['stem_completely_paired'],
			entry['stem_paired_with_itself'],
			entry['stem_paired_nucleotides'],
			entry['stem_unpaired_nucleotides'],
			entry['loop_completely_unpaired'],
			entry['loop_longest_unpaired_stretch'],
			entry['loop_paired_nucleotides'],
			entry['loop_unpaired_nucleotides'],
			entry['linker_completely_unpaired'],
			entry['linker_completely_unpaired'],
			entry['linker_unpaired_nucleotides'],
			entry['mfe'],
			entry['position_on_origin'],
			entry['structure_of_origin_position'],
			entry['linker_name'],
			entry['stem_name']
		)

	# Write output string to file
	with open(output, 'w') as file:
		file.write(output_string)


#######################
#######################
###    functions    ###
#######################
#######################
# Print error message and exit script	
def errx(message):
    print(message)
    exit(1)
    
def has_repeated_char(sequence,char,length):
	pattern = re.compile(f'{re.escape(char)}{{{length}}}')
	return bool(pattern.search(sequence))

def has_number_char(sequence,char,amount):
	hits = sequence.count(char)
	return True if hits >= amount else False

def get_longest_consecutive_char_chain(sequence,char):
	max 			= 0
	current_length 	= 0

	for current_char in sequence:
		if current_char == char:
			current_length += 1
			if current_length > max:
				max = current_length
		else:
			current_length = 0

	return(max)

##########################
### starts main script ###
##########################
main(sequence_file=args.sequence_file,
	output=args.output,
	length_loop=args.length_loop)