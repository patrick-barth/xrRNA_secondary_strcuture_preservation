#!/usr/bin/env python

import random
import os
from Bio import Seq
from Bio import SeqIO
from Bio.SeqUtils import nt_search

##############################################
# Generates a chosen amount of DNA sequences #
#  which do not contain specific motifs and  #
#  have a distinct GC-content                #
##############################################

__author__ = "Patrick Barth" 
__email__ = "Patrick.Barth@biologie.uni-regensburg.de"

#######################
###    variables    ###
#######################

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


def generate_random_seq(seq_len = 100, gc_content = 50, illegal_seqs = None, use_default_illegal_seqs = None, illegal_seq_file = None):
	remove_seqs = use_default_illegal_seqs or illegal_seq_file
	seq_len = int(seq_len)
	gc_content = float(gc_content)
	collected_illegal_seqs = []
	
	# Check if a file of illegal seqs is provided
	if illegal_seq_file:
		# parse file containing illegal sequences. Maybe check if it's a fasta file or a seq file
		collected_illegal_seqs.extend(illegal_seq_file)
		errx('not yet implemented')
	# add default illegal seqs when demanded
	if use_default_illegal_seqs:
		collected_illegal_seqs.extend(base_illegal_sequences)

	if illegal_seqs:
		collected_illegal_seqs.extend(illegal_seqs)
	
	# If some seqs shall not appear in returned sequence then they are transformed to a Seq object
	if remove_seqs:
		for i in range(len(collected_illegal_seqs)):
			collected_illegal_seqs[i] = Seq.Seq(collected_illegal_seqs[i])

	current_seq = Seq.Seq('')

	while len(current_seq) <= seq_len:
		## //TODO: Look at the previous distribution and slightly change the modifier:
		### See how the probability needs to change in order to get to the same GC ratio after a certain motif was removed
		current_seq += get_nucleotide(gc_content,None)

		# if the wanted length is reached it is checked if any unwanted sequences are present (sense and reverse complementary)
		if remove_seqs and len(current_seq) >= seq_len:
			for forward_seq in collected_illegal_seqs:
				reverse_seq = forward_seq.reverse_complement()
				for seq in [forward_seq,reverse_seq]:						

					motif_pos = nt_search(str(current_seq),str(seq))
					while not len(motif_pos) <= 1: 
						current_seq = current_seq[:motif_pos[1]] + current_seq[motif_pos[1] + len(seq):]
						motif_pos = nt_search(str(current_seq),str(seq))

	return(current_seq)

#######################
###    functions    ###
#######################
    
def errx(message):
    print(message)
    exit(1)

def get_nucleotide(gc_chance,omitted_nuc):
	gc_chance = gc_chance/100
	rnd_number  = random.uniform(0,1)
	rnd_number2 = random.uniform(0,1)
	if rnd_number < gc_chance:
		if rnd_number2 < 0.5:
			return 'C'
		else:
			return 'G'
	else:
		if rnd_number2 < 0.5:
			return 'A'
		else:
			return 'T'