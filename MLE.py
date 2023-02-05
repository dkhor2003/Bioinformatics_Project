import math

# open fasta file for extracting 10 nuleotide positions
with open("MA0001.1.fa", "r") as seq:

		#Initialize an empty list for storing all 10 nucleotides sequences
        seq_list =[]
        for lines in seq:

				# If line does not start with ">", it is a sequence
                if not lines.startswith(">"):
                        for character in lines:

								# If ASCII value of character is between 65 and 90, it is a capital letter
                                if ord(character) in range(65, 91):

										# Find index of first capital letter in sequence then break out of loop
                                        seq_index = lines.find(character)
                                        break

						# Extract the 10 capital letter nucleotides
                        sequence = lines[seq_index : seq_index + 10]

						# Append the 10 nucleotide sequence into a list
                        seq_list.append(sequence)

# Initialize four dictionaries for storing the counts of each type of nucleotides at positions 1 - 10
NA = {"1" : 0, "2" : 0, "3" : 0, "4" : 0, "5" : 0, "6" : 0, "7" : 0, "8" : 0, "9" : 0, "10" : 0}
NC = {"1" : 0, "2" : 0, "3" : 0, "4" : 0, "5" : 0, "6" : 0, "7" : 0, "8" : 0, "9" : 0, "10" : 0}
NG = {"1" : 0, "2" : 0, "3" : 0, "4" : 0, "5" : 0, "6" : 0, "7" : 0, "8" : 0, "9" : 0, "10" : 0}
NT = {"1" : 0, "2" : 0, "3" : 0, "4" : 0, "5" : 0, "6" : 0, "7" : 0, "8" : 0, "9" : 0, "10" : 0}

# Counts the number of each type of nucleotide present at positions 1 - 10 of each sequences stored in the list
for elements in seq_list:
        for nuc in range(0, len(elements)):
                if elements[nuc] == 'A':
                        NA[str(nuc + 1)] = NA.get(str(nuc + 1)) + 1
                if elements[nuc] == 'C':
                        NC[str(nuc + 1)] = NC.get(str(nuc + 1)) + 1
                if elements[nuc] == 'G':
                        NG[str(nuc + 1)] = NG.get(str(nuc + 1)) + 1
                if elements[nuc] == 'T':
                        NT[str(nuc + 1)] = NT.get(str(nuc + 1)) + 1


# Function that computes the MLE estimator for the probability of A, C, T, and G for positions 1 - 10
def MLE_function(Na, Nc, Ng, Nt):
        MLE = {}
        MLE["Pa"] = "%.3f"%(Na / (Na + Nc + Ng + Nt))
        MLE["Pc"] = "%.3f"%(Nc / (Na + Nc + Ng + Nt))
        MLE["Pg"] = "%.3f"%(Ng / (Na + Nc + Ng + Nt))
        MLE["Pt"] = "%.3f"%(Nt / (Na + Nc + Ng + Nt))
        return MLE


# Function that computes the sequence logo values for each nucleotides for positions 1 - 10
def Sequence_Logo(Na, Nc, Ng, Nt):
        prob = MLE_function(Na, Nc, Ng, Nt)
        SL = {}
        r = Ri(Na, Nc, Ng, Nt)
        SL["A"] = "%.3f"%(float(prob.get("Pa")) * r)
        SL["C"] = "%.3f"%(float(prob.get("Pc")) * r)
        SL["G"] = "%.3f"%(float(prob.get("Pg")) * r)
        SL["T"] = "%.3f"%(float(prob.get("Pt")) * r)
        
        return SL


# Computes Ri, which is needed to calculate the sequence logo values. 
def Ri(Na, Nc, Ng, Nt):
        SE = Shannon_Entropy(Na, Nc, Ng, Nt)
        r = math.log2(4) - SE + (3 / (2 * 97 * math.log(2)))
        return r



# Computes the Shannon entropy. In cases of which the MLE estimator (x) is 0, the limit as x--> 0 of x log(x) is 0.
def Shannon_Entropy(Na, Nc, Ng, Nt):
        MLE_prob = MLE_function(Na, Nc, Ng, Nt)
        sum = 0
        for MLEs in MLE_prob:
                if not float(MLE_prob.get(MLEs)) ==  0.000:
                        sum = sum + (float(MLE_prob.get(MLEs)) * math.log2(float(MLE_prob.get(MLEs))))
                
        return -sum


# Prints out the 10 by 4 matrix of the MLE estimators for the probability of each nucleotide [A, C, G, T] at positions
# 1 - 10. The format is as follows, where P1A indicates the MLE estimator of the probability of nucleotide A at position 1
#  
#   P1A    P1C    P1G    P1T
#   P2A    P2C    P2G    P2T
#   P3A    P3C    P3G    P3T
#   P4A    P4C    P4G    P4T
#   P5A    P5C    P5G    P5T
#   P6A    P6C    P6G    P6T
#   P7A    P7C    P7G    P7T
#   P8A    P8C    P8G    P8T
#   P9A    P9C    P9G    P9T
#   P10A   P10C   P10G   P10T
for i in range(1, 11):
        MLE_by_position = MLE_function(NA.get(str(i)), NC.get(str(i)), NG.get(str(i)), NT.get(str(i)))
        print(MLE_by_position["Pa"] + "\t" + MLE_by_position["Pc"] + "\t" + MLE_by_position["Pg"] + "\t" + MLE_by_position["Pt"])

print("\n")

# Prints out the 10 by 4 matrix of the sequence logo values for each nucleotide [A, C, G, T] 
# at positions 1 - 10
# The format is as follows, where 1A indicates the sequence logo value of nucleotide A at position 1
#  
#   1A    1C    1G    1T
#   2A    2C    2G    2T
#   3A    3C    3G    3T
#   4A    4C    4G    4T
#   5A    5C    5G    5T
#   6A    6C    6G    6T
#   7A    7C    7G    7T
#   8A    8C    8G    8T
#   9A    9C    9G    9T
#   10A   10C   10G   10T
for j in range(1, 11):
        seq_logo = Sequence_Logo(NA.get(str(j)), NC.get(str(j)), NG.get(str(j)), NT.get(str(j)))
        print(seq_logo["A"] + "\t" + seq_logo["C"] + "\t" + seq_logo["G"] + "\t" + seq_logo["T"])
