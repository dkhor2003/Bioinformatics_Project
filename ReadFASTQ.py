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
# [KSD] You should give your functions names that are more self-explanatory.
def common_function(UMI_line, UMI_dictionary):
    if "ATTGCGCAATG" in UMI_line: # [KSD] This 11bp tag must be at the 5' end of the read.

        # If there are 8 nucleotides after the 11bp tag ATTGCGCAATG
        if UMI_line.index("ATTGCGCAATG") + 11 + 8 <= len(UMI_line):

            # Increments the "count" key's value in the dictionary. This corresponds to the total count of all
            # 8bp UMIs found in the file
            UMI_dictionary["count"] = UMI_dictionary.get("count") + 1

            # Gets the index of the first nucleotide letter in the 8bp UMI
            UMI_index = UMI_line.index("ATTGCGCAATG") + 11

            # Extract 8bp UMI from the line
            UMI = UMI_line[UMI_index: UMI_index + 8]

            # If the 8bp UMI is a new sequence that is not present as one of the keys in the dictionary,
            # add a new key value pair to the dictionary with the sequence of 8bp UMI as the key and 1 as its value
            if UMI not in UMI_dictionary:
                UMI_dictionary[UMI] = 1

            # The 8bp UMI already has a key with the same sequence in the dictionary, increment its value by 1
            else:
                UMI_dictionary[UMI] = UMI_dictionary.get(UMI) + 1


# Opens "combined.read1.fastq" for reading
with open("../combined.read1.ss.fastq", "r") as read_file:	# [KSD] edited

    # Initialize a dictionary with a key called "count" to store the total count of all 8bp UMIs in the file
    # as its value, which is initialized as 0 at the start
    UMI_dict = {"count": 0}

    # Repeat the file reading process for 67061 times, which is the floored number I get from dividing the total
    # number of lines in one of the files by 40000.
    #for z in range(67061):			# [KSD] removed; this does not help unless it is parallelized
    #    fileReader(read_file, UMI_dict)

    # By setting less to True, I am reading the remaining lines in the files, which is less than 40000 lines.
    fileReader(read_file, UMI_dict, less=True)

    # Open a text file for writing results
    with open("UMI.txt", "w") as wFile:
        for seq in UMI_dict:
            if seq != "count":
                wFile.write(seq + "\t" + str(UMI_dict.get(seq)) + "\n")
        wFile.write("Total count: " + str(UMI_dict.get("count")) + "\n")
