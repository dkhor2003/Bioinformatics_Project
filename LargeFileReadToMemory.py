# Function that reads the consequent 40000 lines of the given files if 'less' equals False;
# This function will read the remaining lines in the given files if 'less' equals True
def file_reader(file1, file2, barcode_dictionary, less=False):
    if not less:
        for i in range(1, 40001):
            index1 = file1.readline().rstrip()
            index2 = file2.readline().rstrip()
            # If i % 4 == 2, the current lines consists of the index
            if i % 4 == 2:
                # Increment the "reads" key's value in the dictionary by 1. This corresponds to the total count of
                # barcodes that can be joined, which also indicates the number of indexes in each file.
                barcode_dictionary["reads"] = barcode_dictionary.get("reads") + 1
                common_function(index1, index2, barcode_dictionary)
    else:
        j = 0
        for lines1, lines2 in zip(file1, file2):
            j = j + 1
            # If j % 4 == 2, the current lines consists of the index
            if j % 4 == 2:
                # Increment the "reads" key's value in the dictionary by 1. This corresponds to the total count of
                # barcodes that can be joined, which also indicates the number of indexes in each file.
                barcode_dictionary["reads"] = barcode_dictionary.get("reads") + 1
                index1 = lines1.rstrip()
                index2 = lines2.rstrip()
                common_function(index1, index2, barcode_dictionary)


# Function that joins the index1 and index2 together into a barcode, then search for it in the given dictionary
def common_function(index1, index2, barcode_dictionary):
    barcode = index1 + index2
    # if the joined barcode is in the dictionary, increment its count by 1
    if barcode in barcode_dictionary:
        barcode_dictionary[barcode] = barcode_dictionary.get(barcode) + 1
    # if the joined barcode contains sequence "GGGGGGGG", increment the "GGGGGGGG" key's value in the dictionary by 1
    if "GGGGGGGG" in barcode:
        barcode_dict["GGGGGGGG"] = barcode_dict.get("GGGGGGGG") + 1
    # if the joined barcode contains the nucleotide "N", increment the "N" key's value in the dictionary by 1
    if "N" in barcode:
        barcode_dict["N"] = barcode_dict.get("N") + 1


# Open "Smartseq3.annotation.txt" to read and extract all the unique barcodes, then put them in a dictionary as keys
# with initial values of 0. I also include a key called "reads", a key called "GGGGGGGG", and a key called "N"
# for checking purposes.
with open("../Smartseq3.annotation.txt", "r") as barcode_list:	# [KSD] edited
    barcode_dict = {}
    for lines in barcode_list:
        if lines.startswith("index1"):	# [KSD] For efficiency, discard the first line to avoid check for EVERY line.
            continue
        seq_info = lines.split()
        BC = seq_info[2]
        barcode_dict[BC] = 0
    barcode_dict["reads"] = 0
    barcode_dict["GGGGGGGG"] = 0
    barcode_dict["N"] = 0

    # Open both "combined.index1.fastq" and "combined.index2.fastq" for reading to check if the concatenated barcode
    # exists in the initialized dictionary above.
    with open("../combined.index1.ss.fastq", "r") as index1_file, open("../combined.index2.ss.fastq", "r") as index2_file:	# [KSD] edit
        # Repeat the file reading process for 67061 times, which is the floored number I get from dividing the total
        # number of lines in one of the files by 40000.
        #for loop in range(67061):	# [KSD] edited (it is not a good idea to hard-code numbers that depend on eternal info, like file size)
        #    file_reader(index1_file, index2_file, barcode_dict)

        # By setting less to True, I am reading the remaining lines in the files, which is less than 40000 lines.
        file_reader(index1_file, index2_file, barcode_dict, less=True)	# [KSD] Read the whole file this way

        # A variable for counting the total number of unique barcodes in the file as presented in the file
        # "Smartseq3.annotation.txt"
        unique_read_counts = 0

        # Open a text file for writing results
        with open("barcode_reads.txt", "w") as barcode_reads:
            # Write all the keys (unique barcode) except "reads" and "GGGGGGGG" in the dictionary into the writing
            # file as well as their respective counts.
            for elements in barcode_dict:
                if elements != "reads" and elements != "GGGGGGGG" and elements != "N":
                    barcode_reads.write(elements + "\t" + str(barcode_dict.get(elements)) + "\n")
                    unique_read_counts = unique_read_counts + barcode_dict.get(elements)

            barcode_reads.write("Total reads: " + str(barcode_dict.get("reads")) + "\n")
            barcode_reads.write("Total unique read counts: " + str(unique_read_counts) + "\n")	# [KSD] Had to read the code to understand what this was.
            barcode_reads.write("Reads containing GGGGGGGG: " + str(barcode_dict.get("GGGGGGGG")) + "\n")
            barcode_reads.write("Reads containing N: " + str(barcode_dict.get("N")) + "\n")
