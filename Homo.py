# Opens "alignments.txt" file to read
with open("alignments.txt", "r") as align:
	
	# Initialize a dictionary containing all the 25 possible combinations of alignments as the keys, with values
	# of 0. "AC" means nucleotide A in sequence 1 is aligned with nucleotide C in sequence 2. "-N" means a gap in sequence 
	# 1 is aligned with N, where N = {A,C,G,T}. We consider -A, -C, -G, and -T all as -N because the alignment does not 
	# differentiate between what nucleotide the gap is aligned to. The same applieas to N-.
	align_dict = {"AA" : 0, "AC" : 0, "AG" : 0, "AT" : 0, 
	              "CA" : 0, "CC" : 0, "CG" : 0, "CT" : 0, 
                  "GA" : 0, "GC" : 0, "GG" : 0, "GT" : 0, 
                  "TA" : 0, "TC" : 0, "TG" : 0, "TT" : 0, 
                  "-N" : 0, "N-" : 0, "--" : 0}
	
	# Loop through all the alignments to count the 19 possible combinations of alignments and increment the corresponding 
	# key's value in the dictionary. 
	for z in range(10):
		
		# Get sequence 1 of alignment
		seq1 = align.readline().rstrip()

		# Get sequence 2 of alignment
		seq2 = align.readline().rstrip()
		startIndex = 0

		# If sequence 1 starts with '-', it is an end gap
		if seq1.startswith("-"):
			
			# Loop through the sequence and find the index after the end gap ends. Set this this to be the start index. We 
			# don't want to include alignments with end gaps in our counting since we expect it to be there in semi-global 
			# alignment
			for i in range(len(seq1)):
				if seq1[i] != "-":
					startIndex = i
					break

		# If sequence 2 starts with '-', it is an end gap
		if seq2.startswith("-"):
			
			# Loop through the sequence and find the index after the end gap ends. Set this this to be the start index. We 
			# don't want to include alignments with end gaps in our counting since we expect it to be there in semi-global 
			# alignment
			for j in range(len(seq2)):
				if seq2[j] != "-":
					startIndex = j
					break

		# Find the end gap starting index in the other sequence. For example, if sequence 1 has an end gap at the start, 
		# then sequence 2 will have an end gap at the end. 
		endIndex = max(seq1.find("---"), seq2.find("---"))

		# Loop through each nucleotide or gap between the start index and end index of both sequence 1 and 2, which excludes 
		# the end gaps.Counts the corresponding aligned characters between sequence 1 and 2 at a specific position and 
		# increment the corresponding key's value in the dictionary by 1. 
		for i in range(startIndex, endIndex):
			if seq2[i] == ".":
				combined = 2 * seq1[i]
			elif seq1[i] == "-" and seq2[i] != "-": # Gap in sequence 1 aligned to a nucleotide in sequence 2
				combined = "-N"
			elif seq1[i] != "-" and seq2[i] == "-": # Gap in sequence 2 aligned to a nucleotide in sequence 1
				combined = "N-"
			else:
				combined = seq1[i] + seq2[i]
			align_dict[combined] = int(align_dict.get(combined)) + 1


# Returns a dictionary of the MLEs of the probabilities of each possible alignments using the MLE equation derived based on 
# a homologous model 
def Homo_MLE_function(dictionary):
	MLE = {}
	total = 0
	for elements in dictionary:
		total = total + dictionary.get(elements)
	for elements in dictionary:
		
		# If the key is -N or N-, we need to divide the MLE by 4 because it accounts for four different nucleotide 
		# A, C, G, and T
		MLE[elements] = dictionary.get(elements) / total
	return MLE


homo = Homo_MLE_function(align_dict)

# Prints the 5x5 matrix of the MLE probabilities of each possible alignment in the format:
#
#	PAA  PAC  PAG  PAT  PA-
#	PCA  PCC  PCG  PCT  PC-
#	PGA  PGC  PGG  PGT  PG-
#	PTA  PTC  PTG  PTT  PT-
#	P-A  P-C  P-G  P-T  P--
#
#	where P-A,  P-C,  P-G, and P-T will all have the same MLE as each other since the alignment does not differentiate 
#	between whether the gap in sequence 1 will be aligned to what type of nucleotide. 
#	The same applies to PA-, PC-, PG-, and PT-. 
if __name__ == '__main__':
	for i in ["A", "C", "G", "T", "-"]:
		for j in ["A", "C", "G", "T", "-"]:
			if i == "-" and j != "-":
				combined = "-N"
				
				# Although not correct, the MLE for P-A, P-C, P-G, and P-T are each the result of taking P-N divided by 4
				# This is done so to make sure the sum of all MLEs of the 5 by 5 matrix equals to 1
				print("%.3f"%(homo.get(combined) / 4), end = "\t")
			elif i != "-" and j == "-":
				combined = "N-"
				
				# Although not correct, the MLE for PA-, PC-, PG-, and PT- are each the result of taking PN- divided by 4
				# This is done so to make sure the sum of all MLEs of the 5 by 5 matrix equals to 1
				print("%.3f"%(homo.get(combined) / 4), end = "\t")
			else:
				combined = i + j
				print("%.3f"%(homo.get(combined)), end = "\t")
		print()

