#!/usr/bin/python

import numpy as np
import math
from Bio.Align import substitution_matrices
from scipy import optimize

q_aa = {'A' : 0.074, 'Q' : 0.034, 'L' : 0.099, 'S' : 0.057,
	'R' : 0.052, 'E' : 0.054, 'K' : 0.058, 'T' : 0.051,
	'N' : 0.045, 'G' : 0.074, 'M' : 0.025, 'W' : 0.013,
	'D' : 0.053, 'H' : 0.026, 'F' : 0.047, 'Y' : 0.033,
	'C' : 0.025, 'I' : 0.068, 'P' : 0.039, 'V' : 0.073}
blosum62 = substitution_matrices.load('BLOSUM62')

# Cutoff value
small_value = 1e-5
##
# Function for which we seak zeros.
#
# @param x	current value of lambda
# @param q	probabilities of each letter in alphabet
# @param scores	scores of all possible mismatches and matches in double dictionary
def func(x, q = q_aa, scores = blosum62):
	sum = 0
	for aa1 in q.keys():
		for aa2 in q.keys():
			sum += q.get(aa1) * q.get(aa2) * math.exp(x * scores[aa1][aa2])# a function of q[aa1], q[aa2], scores[aa1][aa2], and x (lambda)
	return(sum - 1)


# Function that gets the non_zero root of a function given a cutoff value (minimum_cutoff) by incorporating the bisection 
# method of scipy.optimize. Here, 'start' and 'end' is the interval of x to search, and 'skip' is the value of how 
# much to skip between each iteration of x. 
def find_non_zero_root(minimum_cutoff, function, start, end, skip):
	
	# Initialize two variables to None at start, where "positive" indicates a x value that results in the function 
	# being positive, and "negative" indicates a x value that results in the function being negative.
	positive = negative = None
	
	# Loop the value of x from 'start' to 'end' with 'skip' intervals
	for i in np.arange(start, end, skip):
		if negative != None and positive != None:
			
			# Get the root of f(x) for x between 'negative' and 'positive'
			res = optimize.bisect(function, a = negative, b = positive)
			
			# If the magnitude of root is greater than the specified cutoff value, return the root
			if abs(res) > minimum_cutoff:
				return res
		
		# If f(x) is negative, set 'negative' as x
		if function(i) < 0:
			negative = i
			
		# If f(x) is positive, set 'positive' as x
		if function(i) > 0:
			positive = i


# Get the non_zero_lambda using the desired cutoff value. 
res = find_non_zero_root(small_value, func, -5, 5, 0.1)
if __name__ == '__main__':
	print(res)


