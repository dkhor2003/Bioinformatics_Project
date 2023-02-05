from numpy.random import default_rng


# Function that randomly generates the same number of nucleotide sequences as the "count" key's value in the given
# dictionary, each with a length specified in the parameter called seq_len. The randomly generated sequences and their
# corresponding count are returned in the form of a dictionary
# [KSD] It would be faster to simulate UMI counts as draw from Multinomial distribution with all probabilities equal to 1/4**8.
def random_data_simulator(UMI_dictionary, seq_len):
    # Construct a generator
    rng = default_rng()

    # Initialize a dictionary with a key called "count" to store the total count of UMIs generated
    random_dict = {"count": 0}

    # Randomly generates the same number of nucleotide sequences as the "count" key's value in the given
    # dictionary, each with a length specified in the parameter called seq_len. Each of the four nucleotide
    # letter has the same probability being chosen to put into the sequence.
    for loop in range(UMI_dictionary.get("count")):
        rand = "".join(rng.choice(['A', 'C', 'G', 'T'], size=seq_len, p=[0.25, 0.25, 0.25, 0.25]))
        random_dict["count"] = random_dict.get("count") + 1
        if rand not in random_dict:
            random_dict[rand] = 1
        else:
            random_dict[rand] = random_dict.get(rand) + 1

    return random_dict


# Computes the test statistic of a given dataset that is the form of a dictionary and returns the value
def getTestStat(prob_dict):
    test_value = 0
    for probs in prob_dict:
        if probs != "count":
            test_value = test_value + (int(prob_dict.get(probs)) * float(prob_dict.get(probs) / prob_dict.get("count")))
    # [KSD] More computationally efficient to compute sum of squared counts and divide by count^2.
    # [KSD] On the other hand, such calculation is prone to numerical overflow.
    return test_value / prob_dict.get("count")


# Computes the p-value by comparing the multiple simulated test statistics with the observed test statistic
def get_pValue(UMI_dictionary, simulation_times):
    # Computes the test statistic of the observed data
    data_test_stat = getTestStat(UMI_dictionary)
    reference_list = []

    # Generates a list containing all different simulated test statistics with a size of simulation_times
    for i in range(simulation_times):
        simulation_dict = random_data_simulator(UMI_dictionary, 8)
        simulation_test_stat = getTestStat(simulation_dict)
        reference_list.append(simulation_test_stat)
    c = 0

    # Loop through all simulated test statistics in the list and compare it with the observed test statistic. If
    # the simulated test statistic is greater than the observed test statistic, the counter variable is incremented by 1
    for value in reference_list:
        if value > data_test_stat:
            c += 1

    # If c = 0, this means that the p-value is 0. However, p-values generally can't be 0 because this will mean
    # that you are 100% sure that there can never be a data with such a value as you observed. Thus, to prevent from
    # getting a p-value of 0, the numerator and denominator is incremented by 1 in case if calculated p-value is 0
    if (c / simulation_times) == 0:
        return (c + 1) / (simulation_times + 1)
    else:
        return c / simulation_times


# Function that reads the consequent 40000 lines of the given file if 'less' equals False;
# This function will read the remaining lines in the given file if 'less' equals True
def fileReader(fileName, UMIdict, less=False):
    if not less:
        for i in range(1, 40001):
            UMI_line = fileName.readline().rstrip()

            # If i % 4 == 2, the current lines consists of the read sequence
            if i % 4 == 2:
                common_function(UMI_line, UMIdict)
    else:
        j = 0
        for lines in fileName:
            j = j + 1

            # If j % 4 == 2, the current lines consists of the read sequence
            if j % 4 == 2:
                common_function(lines, UMIdict)


# Function that checks if the 11bp tag ATTGCGCAATG is in the given line and modifies the dictionary accordingly
def common_function(UMI_line, UMI_dictionary):
    if "ATTGCGCAATG" in UMI_line:    # [KSD] This is a mistake. You need the 11bp tag at the 5' end of the read.

        # Checks if there are 8 nucleotides after the 11bp tag ATTGCGCAATG
        if UMI_line.index("ATTGCGCAATG") + 11 + 8 <= len(UMI_line):   # [KSD] Good check for possibly running of line end.

            # Increments the "count" key's value in the dictionary. This corresponds to the total count of all
            # 8bp UMIs found in the file
            UMI_dictionary["count"] = UMI_dictionary.get("count") + 1

            # Gets the index of the first nucleotide letter in the 8bp UMI
            UMI_index = UMI_line.index("ATTGCGCAATG") + 11

            # Extract 8bp UMI from the line
            UMI = UMI_line[UMI_index: UMI_index + 8]

            # If UMI does not contain any ambiguous nucleotide 'N'，we extract it
            if 'N' not in UMI:

                # If the 8bp UMI is a new sequence that is not present as one of the keys in the dictionary,
                # add a new key value pair to the dictionary with the sequence of 8bp UMI as the key and 1 as its value
                if UMI not in UMI_dictionary:
                    UMI_dictionary[UMI] = 1

                # The 8bp UMI already has a key with the same sequence in the dictionary, increment its value by 1
                else:
                    UMI_dictionary[UMI] = UMI_dictionary.get(UMI) + 1


# Opens "combined.read1.fastq" for reading
with open("../combined.read1.ss.fastq", "r") as read_file: # [KSD] edited
    # Initialize a dictionary with a key called "count" to store the total count of all 8bp UMIs in the file
    # as its value, which is initialized as 0 at the start
    UMI_dict = {"count": 0}

    # Repeat the file reading process for 10 times to sample first 100000 reads in the file
    #for z in range(10):             # [KSD] edited
    fileReader(read_file, UMI_dict)  # [KSD] edited

    # Computes the p-value by inputting the observed data and simulate a random dataset with same counts as the observed
    # dataset for 300 times
    pValue = get_pValue(UMI_dict, 300)

    # Open a text file for writing results
    with open("pValue.txt", "w") as wFile:
        for seq in UMI_dict:
            if seq != "count":
                wFile.write(seq + "\t" + str(UMI_dict.get(seq)) + "\n")
        wFile.write("Total count: " + str(UMI_dict.get("count")) + "\n")
        wFile.write("P-value: " + str(pValue) + "\n")
