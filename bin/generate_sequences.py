#!/usr/bin/env python

#################################################
# Generates sequences from different modules 	#
#  with potentially random linker sequences		#
#################################################

import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from random_sequence_generator import generate_random_seq
from Bio.SeqUtils import nt_search

parser = argparse.ArgumentParser()
parser.add_argument('--xrRNA_file', 	'-x', type=str,	help='Input file with xrRNA motifs')
parser.add_argument('--spacer_file', 	'-s', type=str, default=None,	help='Input file with spacer sequences')
parser.add_argument('--IRES_file', 	'-i', type=str,	help='Input file with IRES sequences')
parser.add_argument('--CDS_file', 	'-c', type=str,	help='Input file with CDS sequences')
parser.add_argument('--only_random_spacers', 	'-r', action=argparse.BooleanOptionalAction, default=False,	help='If true only random spacers will be used')
parser.add_argument('--only_specific_spacers', 	'-a', action=argparse.BooleanOptionalAction, default=False,	help='If true only spacers from delivered file will be used')
parser.add_argument('--percent_specific_spacers', 	'-p', type=float, default=50,	help='Percentage of spacers to be from file. Only works when both - specific and random - spacers are to be used')
parser.add_argument('--spacer_min_len', 	'-b', type=int, default=10,	help='Minimum length of spacer between xrRNA and IRES')
parser.add_argument('--spacer_max_len', 	'-d', type=int, default=1000,	help='Maximum length of spacer between xrRNA and IRES')
parser.add_argument('--spacer_fixed_length', 	'-f', type=int, default=None,	help='If value is provided the spacer will have a fixed length')
parser.add_argument('--number_sequences', 	'-n', type=int, default=100,	help='Number of sequences to be generated')
parser.add_argument('--output', 	'-o', type=str,	help='Output file')
args = parser.parse_args()

########################
########################
###    parameters    ###
########################
########################

base_illegal_sequences = [
		"ATTTAAAT",
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

spacer_5_prime_end = 'GCAA'
spacer_3_prime_end = 'ACTA'

##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main(xrRNA_file,spacer_file,IRES_file,CDS_file,only_random_spacers,only_specific_spacers,percent_specific_spacers,spacer_min_len,spacer_max_len,spacer_fixed_length,number_sequences,output):
	#TODO: Check if there is a more elegant way to ensure that only one can be used
	if only_random_spacers and only_specific_spacers:
		errx('Error: Both options - only_random_spacers and only_specific_spacers - have been set to True.\nThese cannot be used at the same time as they exclude each other.')
	if only_specific_spacers and spacer_file is None:
		errx('Error: You selected that only random spacers are to be used but delivered no file with spacer sequences.')
	if percent_specific_spacers < 0 or percent_specific_spacers > 100:
		errx('Percentage of specific spacers ot be used is not within allowed values.\n--percent_specific_spacers needs to be set between 0 and 100')
	
	# Read fasta files
	xrRNA_sequences 	= import_fasta(xrRNA_file)
	if not only_random_spacers or spacer_file is not None:
		spacer_sequences 	= import_fasta(spacer_file)
	ires_sequences 		= import_fasta(IRES_file)
	cds_sequences 		= import_fasta(CDS_file)

	# Initiate Sequence collection
	collected_sequences = []

	i = 0
	while i < number_sequences:

		tmp_xrRNA_seq = random.choice(xrRNA_sequences)
		# Check if spacer sequences have been read
		if 'spacer_sequences' in locals() and only_specific_spacers:
			tmp_spacer_entry 	= random.choice(spacer_sequences)
			tmp_spacer_seq 		= tmp_spacer_entry.seq
			tmp_spacer_id 		= tmp_spacer_entry.id
		else:
			# Get fixed spacer length if desired. Else random length is chosen
			spacer_seq_length = spacer_fixed_length if spacer_fixed_length else random.randint(spacer_min_len,spacer_max_len)
			
			if 'spacer_sequences' in locals() and not only_random_spacers:
				if get_random_chance(percent_specific_spacers):
					tmp_spacer_entry 	= random.choice(spacer_sequences)
					tmp_spacer_seq 		= tmp_spacer_entry.seq
					tmp_spacer_id 		= tmp_spacer_entry.id
				else:
					tmp_spacer_seq 	= generate_random_seq(seq_len = spacer_seq_length, illegal_seqs = base_illegal_sequences)
					tmp_spacer_id 	= 'random'
			else:
				tmp_spacer_seq 	= generate_random_seq(seq_len = spacer_seq_length, illegal_seqs = base_illegal_sequences)
				tmp_spacer_id 	= 'random'
		# Add adapter necessary for cloning
		tmp_spacer_seq = Seq(spacer_5_prime_end) + tmp_spacer_seq + Seq(spacer_3_prime_end)
		tmp_ires_seq 	= random.choice(ires_sequences)
		tmp_cds_seq 	= random.choice(cds_sequences)

		## assemble SeqRecord object
		# Get positions for every module
		xrRNA_start_pos = 0
		xrRNA_end_pos = len(tmp_xrRNA_seq.seq)
		spacer_start_pos = len(tmp_xrRNA_seq.seq) + 1
		spacer_end_pos = spacer_start_pos + len(tmp_spacer_seq)
		ires_start_pos = spacer_end_pos + 1
		ires_end_pos = ires_start_pos + len(tmp_ires_seq)
		cds_start_pos = ires_end_pos +1
		cds_end_pos = cds_start_pos + len(tmp_cds_seq)
		# TODO: possibility to add spacer between IRES and CDS
		id = str(i) + '_' +tmp_xrRNA_seq.id + ';' + tmp_spacer_id + ';' + tmp_ires_seq.id + ';' + tmp_cds_seq.id
		description = (
			f'xrRNA:{xrRNA_start_pos}-{xrRNA_end_pos} '
			f'spacer:{spacer_start_pos}-{spacer_end_pos} '
			f'ires:{ires_start_pos}-{ires_end_pos} '
			f'cds:{cds_start_pos}-{cds_end_pos}'
		)

		sequence = tmp_xrRNA_seq.seq + Seq(tmp_spacer_seq) + tmp_ires_seq.seq + tmp_cds_seq.seq
		sequence = sequence.upper()

		# Check if illegal sequences appear in resulting sequence
		#TODO: find a way to make this work. Current problem is that some of the sequences are found in the fixed modules and
		## that the CDS contains start and termination codons
		#keep_seq = True
		#for forward_illegal_seq in base_illegal_sequences:
		#	forward_illegal_seq = Seq(forward_illegal_seq)
		#	reverse_illegal_sequence = forward_illegal_seq.reverse_complement()
		#	for seq in [forward_illegal_seq,reverse_illegal_sequence]:
		#		if seq in sequence[0:ires_end_pos]:
		#			print(id)
		#			print(description)
		#			print(str(seq) + ' has been found at ' + str(nt_search(str(sequence),str(seq))))
		#			keep_seq = False
		#			break

		#if not keep_seq:
		#	continue

		record = SeqRecord(Seq(sequence),
							id = id,
							description = description)
		
		collected_sequences.append(record)
		i += 1



	with open(output,'w') as output_handle:
		SeqIO.write(collected_sequences, output_handle, "fasta")







#######################
#######################
###    functions    ###
#######################
#######################
# Print error message and exit script	
def errx(message):
    print(message)
    exit(1)
    
def import_fasta(path):
	collected_seqs = []
	for record in SeqIO.parse(path, "fasta"):
		collected_seqs.append(record)
	return collected_seqs

# Provide percent value and receive True based on the percentage
def get_random_chance(percent):
	value = percent/100
	random_value = random.uniform(0,1)
	if random_value < value:
		return True
	else:
		return False

##########################
### starts main script ###
##########################
main(xrRNA_file=args.xrRNA_file,
	spacer_file=args.spacer_file,
	IRES_file=args.IRES_file,
	CDS_file=args.CDS_file,
	only_random_spacers=args.only_random_spacers,
	only_specific_spacers=args.only_specific_spacers,
	percent_specific_spacers=args.percent_specific_spacers,
	spacer_fixed_length=args.spacer_fixed_length,
	spacer_min_len=args.spacer_min_len,
	spacer_max_len=args.spacer_max_len,
	number_sequences=args.number_sequences,
	output=args.output)


