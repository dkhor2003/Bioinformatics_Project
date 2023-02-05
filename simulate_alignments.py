#!/usr/bin/python

from Bio.Align import substitution_matrices
from numpy.random import default_rng
import math

# import find_lambda.py to use variable res
import find_lambda

gap_penalty = -11# choose your gap penalty
estimated_lambda = find_lambda.res # set this from your previous work

num_excursions_to_simulate = 1000
num_pairs_to_simulate = 2000

q_aa = {'A' : 0.074, 'Q' : 0.034, 'L' : 0.099, 'S' : 0.057,
	'R' : 0.052, 'E' : 0.054, 'K' : 0.058, 'T' : 0.051,
	'N' : 0.045, 'G' : 0.074, 'M' : 0.025, 'W' : 0.013,
	'D' : 0.053, 'H' : 0.026, 'F' : 0.047, 'Y' : 0.033,
	'C' : 0.025, 'I' : 0.068, 'P' : 0.039, 'V' : 0.073}

p_g = math.exp(gap_penalty * estimated_lambda)
print("My gap penality is:", p_g)

# probabilities of each alignment pair
aa = q_aa.keys()
aa_probs = [ q_aa[i] * q_aa[j] * (1-2*p_g) for j in aa for i in aa]
aa_del = [ p_g * q_aa[i] for i in aa]
del_aa = [ p_g * q_aa[i] for i in aa]
probs = [*aa_probs, *aa_del, *del_aa]

# the alignment pairs: AA, AC, ..., -A, ..., T-
aa_aa_pairs = [ i + j for i in aa for j in aa ]
aa_del_pairs = [ i + '-' for i in aa ]
del_aa_pairs = [ '-' + j for j in aa ]
pairs = [*aa_aa_pairs, *aa_del_pairs, *del_aa_pairs]

# scoring system
blosum62 = substitution_matrices.load('BLOSUM62')

rng = default_rng()
sim_pairs = rng.choice(pairs, size = num_pairs_to_simulate, p = probs)

R = [0] * 12	# R probabilities
Q = [0] * 12	# Q probabilities

A = 0		# average excursion length
n_excursions = 0
n_pairs = num_pairs_to_simulate
pheight = height = 0
excursion_length = 0

# variable to keep track of highest scoring excursion
max_height = 0
while n_excursions < num_excursions_to_simulate:
	n_pairs -= 1
	sim_pair = sim_pairs[n_pairs]
	excursion_length += 1
	[aa1, aa2] = list(sim_pair)

	if aa1 == '-' or aa2 == '-':
		# keep track of the excursion height
		height += gap_penalty
	else:
		# keep track of the excursion height
		height += blosum62.get((aa1, aa2))
		
		# If the current score is greater than the maximum score, set it as the new maximum score (highest scoring 
		# excursion)
		if height > max_height:
			max_height = height	

	# update all estimates if starting or ending excursion
	if pheight == 0 and height > 0:	# starting
		Q[int(height)] += 1
	if height < 0:			# ending
		R[int(-height)] += 1
		A += excursion_length
		height = 0		# reset
		n_excursions += 1
		excursion_length = 0

	# have completed all required excursions
	if n_excursions == num_excursions_to_simulate:
		break
	pheight = height

	# need to simulate more pairs
	if n_pairs == 0:
		print("Resimulating more data...")
		n_pairs = num_pairs_to_simulate
		sim_pairs = rng.choice(pairs, size = num_pairs_to_simulate, p = probs)
print("Simulated", n_excursions, "excursions")

# output all estimates
sum = 0
for i in range(1,len(Q)):
	Q[i] /= n_excursions
	print("Q[", i, "]:", Q[i])
	sum += Q[i]
print("\\overline R:", 1 - sum)
for i in range(1,len(R)):
	R[i] /= n_excursions
	print("R[", -(i), "]:", R[i])
print("A:", A/n_excursions)

# Computing B based on the equation given in BLAST lecture slides, where: 
#    - prob_Q_negative is the probability that an excursion starts with a negative height
#    - prob_Q_list is the list of probabilities of which the excursion visits k before any other positive value, for 
#      k = 1,..., highest matching score in the given BLOSUM matrix
#    - prob_R_list is the list of probabilities that the excursion ends in state -j, for j = -1,...,gap penalty score
#    - lambda_value is the estimated value of lambda calculated in find_lambda.py
def compute_B(prob_Q_negative, prob_Q_list, prob_R_list, lambda_value):
	t = 1
	for c in range(1, len(prob_R_list)):
		t = t - (prob_R_list[c] * math.exp(-lambda_value * c))

	k = 0
	for d in range(1, len(prob_Q_list)):
		k = k + (d * prob_Q_list[d] * math.exp(d * lambda_value))

	return (prob_Q_negative * t) / ((1 - math.exp(-lambda_value)) * k)

# E-value of highest excursion score: nBe^(-lambda * y), where n is the number of simulated excursions, B is the value 
# gotten from computed_B, lamda is the estimated lambda value calculated in find_lambda.py, and y is the highest excursion 
# height. 
def e_value(B, lambda_value, excursion_number, max_excursion_score):
	return (excursion_number * B * math.exp(-lambda_value * max_excursion_score)) 

# Prints out the highest excursion height
print("Highest excursion height: " + str(max_height))

# Prints out the E-value of the highest excursion score
print("E-value of highest excursion score: " + str(e_value(compute_B(1 - sum,Q,R,estimated_lambda), estimated_lambda,num_excursions_to_simulate,max_height))) 
