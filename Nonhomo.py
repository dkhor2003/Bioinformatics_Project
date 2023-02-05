# Opens "alignments.txt" file to read
with open("alignments.txt", "r") as align:
	
	# Initialize two dictionary containing the 4 nucleotides (A, C, G, T) as their keys with an initial value of 0, which 
	# represents the counts.One for sequence 1 of alignment and another one for sequence 2 of alignment 
	s1 = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}
	s2 = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}

	# Loop through all the alignments to count the appearance of each nucleotide (A, C, G, T) in each sequence of alignment. 
	for z in range(10):

		# Get sequence 1 of alignment
		seq1 = align.readline().rstrip()

		# Get sequence 2 of alignment
		seq2 = align.readline().rstrip()

		# Loop through each nucleotide or gap in both sequence 1 and 2.Counts the each nucleotide that is present in each 
		# sequence and increment the count by 1 correspondingly in each dictionary.s1 for sequence 1 and s2 for sequence 2. 
		# If we encounter a gap ('-'), we ignore it and proceed since gaps don't play a role in nonhomologous model
		for i in range(0, len(seq1)):
			if seq1[i] != "-":
				s1[seq1[i]] = s1.get(seq1[i]) + 1
			if seq2[i] != "-":
				if seq2[i] == ".":
					s2[seq1[i]] = s2.get(seq1[i]) + 1
				else:
					s2[seq2[i]] = s2.get(seq2[i]) + 1

# Returns a dictionary of the MLEs of the probabilities of each possible nucleotide using the MLE equation derived based on 
# a nonhomologous model 
def NonHomo_MLE_function(dict1, dict2):
	total = 0
	MLE ={}
	for element1, element2 in zip(dict1, dict2):
                total = total + dict1.get(element1) + dict2.get(element2)
	MLE["A"] = (dict1.get("A") + dict2.get("A")) / total
	MLE["C"] = (dict1.get("C") + dict2.get("C")) / total
	MLE["G"] = (dict1.get("G") + dict2.get("G")) / total
	MLE["T"] = (dict1.get("T") + dict2.get("T")) / total
	return MLE

nonhomo = NonHomo_MLE_function(s1, s2)

# Prints the MLE probability of each nucleotide in the format:
#
#	PA  PC  PG  PT
#
#	where, for example, PA means the probability of nucleotide A under a nonhomologous model
if __name__ == '__main__':
	for i in nonhomo:
		print("%.3f"%(nonhomo.get(i)), end = "\t")

