#!/bin/python3

from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys

# Import scores.py to use the score dictionary for each possible combination for initializing both the match_dict to 
# be used in globaldd, as well as for getting the gap open penalties and gap extension penalties for both sequence 1
# and sequence 2 of the alignment. 
import scores as sc

f1 = sys.argv[1]		# fasta forward reads
f2 = sys.argv[2]		# fasta reverse reads

# open forward and reverse reads
with open(f1) as fr_file, open(f2) as rr_file:
	
	# Initializing the dictionary to be used for matches and mismatches in globaldd alignment. 
	match_dict = {}
	for i in ["A", "C", "G", "T"]:
		for j in ["A", "C", "G", "T"]:
			match_dict[(i, j)] = sc.scores_dict.get(i + j)
	
	# Variable for the total number of pairs present in the fasta files
	total_pairs = 0
	
	# Variable for the number of mergeable pairs in the fasta files
	mergeable_count = 0
	
	# Variable for the number of readthroughs present in the fasta files
	readthroughs = 0
	
	# taking one pair of reads at a time
	for i, (fr, rr) in enumerate(zip(SeqIO.parse(fr_file, "fastq"), SeqIO.parse(rr_file, "fastq"))):
		total_pairs += 1
		if 'N' in fr.seq or 'N' in rr.seq:
			continue
		# align reads without penalizing end gaps (semi-global alignment)
		alignments = pairwise2.align.globaldd(fr.seq, rr.reverse_complement().seq, 
			#===YOU MUST SET SCORES APPROPRIATELY TO MAKE THIS RUN===,
 			match_dict, sc.scores_dict.get("-N"), sc.scores_dict.get("-N"), sc.scores_dict.get("N-"), 
			sc.scores_dict.get("N-"), penalize_end_gaps = (False, False), one_alignment_only = True)
		
		# If the score of the best alignment is not negative, which means that the pairs are mergeable and have 
		# homology, do the following...
		if not alignments[0][2] < 0:
			
			# Print output in format specified in the homework manual
			print(fr.id + "\t" + rr.id + "\t" + str(alignments[0][2]))
			
			# Increment mergeable pair count variable by 1
			mergeable_count += 1
			
			# If there are end gaps at the start of the first sequence, it means that this is a readthrough. 
			if alignments[0][0].find("---") == 0:
				
				# Increment readthrough count variable by 1
				readthroughs += 1

	print("Total pairs: " + str(total_pairs))
	print("Total mergeable pairs: " + str(mergeable_count))
	print("Total readthroughs: " + str(readthroughs))
