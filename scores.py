import math

# Import both Homo.py and Unhomo.py files to use the MLE dictionaries(Homologous and Nonhomologous) created in both files
import Homo as h1
import Nonhomo as h2

# Function for getting the scores for a particular combination of 'first' and 'second', where 'first' is an element of
# {A, C, G, T, -} in the first sequence, and 'second' is an element of {A, C, G, T, -} in the second sequence. 'homo_dict'
# is the dictionary containing the MLEs for all the possible combinations under the homologous model, and 'nonhomo_dict' is 
# the dictionary containing the MLEs for Pa, Pc, Pg, and Pt under the nonhomologous model. 
def scoring_function(first, second, homo_dict, nonhomo_dict):
	combine = first + second 
	prob1 = float(homo_dict.get(combine))
	
	# If 'first' is a gap and 'second' is a nucleotide (A, C, G, T), its score is just log{PrH(Si1 = -)}
	# If 'second' is a gap and 'first' is a nucleotide (A, C, G, T), its score is just log{PrH(Si2 = -)}
	# 
	# where PrH(Si1 = -) means the probability of the ith site of sequence 1 being a gap under the homology model
	#   and PrH(Si2 = -) means the probability of the ith site of sequence 2 being a gap under the homology model
	if (first == "-" and second != "-") or (first != "-" and second == "-"):
		return math.log(prob1)
	
	# If 'first' and 'second' are both gaps, the score of it will be 0 since the limit of log x as x approaches 0 is 0. 
	elif first == "-" and second == "-":
		return 0
	else:
		prob2 = float(nonhomo_dict.get(first)) * float(nonhomo_dict.get(second))
	return math.log(prob1 / prob2)


scores_dict = {"AA" : 0, "AC" : 0, "AG" : 0, "AT" : 0,
			   "CA" : 0, "CC" : 0, "CG" : 0, "CT" : 0,
			   "GA" : 0, "GC" : 0, "GG" : 0, "GT" : 0,
			   "TA" : 0, "TC" : 0, "TG" : 0, "TT" : 0,
			   "N-" : 0, "-N" : 0}

# Prints out the 5 x 5 matrix of the scores for each possible combination in the following format:
#
# S(A, A)  S(A, C)  S(A, G)  S(A, T)  S(A, -)
# S(C, A)  S(C, C)  S(C, G)  S(C, T)  S(C, -)
# S(G, A)  S(G, C)  S(G, G)  S(G, T)  S(G, -)
# S(T, A)  S(T, C)  S(T, G)  S(T, T)  S(T, -)
# S(-, A)  S(-, C)  S(-, G)  S(-, T)  S(-, -)
#
# where, for example, S(C, G) means the score of nucleotide C in sequence 1 is aligned with nucleotide G in sequence 2
for first in ["A", "C", "G", "T", "-"]:
	for second in ["A", "C", "G", "T", "-"]:
		if first == "-" and second != "-":
			second = 'N'
		elif first != "-" and second == "-":
			first = 'N'
		score = scoring_function(first, second, h1.homo, h2.nonhomo)
		scores_dict[first + second] = score
		if __name__ == '__main__':
			print("%.3f"%(score), end = "\t")
	if __name__ == '__main__':
		print()

