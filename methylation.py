# Import Bio.SeqIO for parsing FASTA file
from Bio import SeqIO

# Function that gets the species name in a read identifier
def get_species(read_identifier):
	l = read_identifier.split("[")
	return  l[len(l) - 1].rstrip("]")
	

# List to hold all the species present in the hexapoda.fsa file 
species = []

# Parse the FASTA file containing hexapoda protein sequences of size ranging from 700 - 1500 amino acids into 
# SeqRecord objects. Then, get the species name and append to the created list with no repeats of the same species name
# in the list. 
for record in SeqIO.parse("hexapoda.fsa", "fasta"):
	s = record.description
	sp = get_species(s)
	if sp not in species:
		species.append(sp)

# Get the total number of hexapoda species found in the FASTA file 
total_species = len(species)

# Function that removes the species name from the list containing species with no evidence of a methylation system and
# adds it into the list containing species that might have a methylation system. 
def add_or_remove_species(species, with_met_sys, without_met_sys):
	if species not in with_met_sys:
		with_met_sys.append(species)
	if species in without_met_sys:
		without_met_sys.remove(species)

# Create a list for storing hexapoda species that might have a methylation system.
species_with_ms = []

# Open the text file containing the BLAST result from blastp with query being the Apis mellifera (honeybee) DNA 
# methyltransferases DNMT1a and DNMT3 RefSeq proteins and the protein database being the database made from the FASTA file 
# hexapoda.fsa. 
with open("hexapoda_blast.txt", "r") as blast:
	
	# Variable for concatenating the current line with the next line in case if the read identifier line does not end 
	# with a "]". 
	# For example: 
	
	#  >XP_026746692.1 LOW QUALITY PROTEIN: DNA (cytosine-5)-methyltransferase PliMCI-like
	#  [Trichoplusia ni]
	#
	#  OR
	#
	#  >XP_045784805.1 DNA (cytosine-5)-methyltransferase PliMCI-like isoform X1 [Maniola
	#  jurtina]
	#
	# In the first case, the species name "[Trichoplusia ni]" is on the second line of the read identifier instead of being
	# in the same line. 
	#
	# In the second case, the species name is separated with the front part of the name in the first line and the latter 
	# part of the name in the second line. 
	#
	# In both cases, we need to concatenate them together to get the full species name
	s_cont = ""
	
	# Boolean variable to keep check whether the read identifier line ends with a "]" or not. True means the line does
	# not end with a "]", False means the line does end with a "]".
	cont = False
	
	# Iterate through each line in the file
	for lines in blast:
		
		# Remove right end whitespaces of each line
		lines = lines.rstrip()
		
		# If there is no need for concatenation
		if cont == False:
			
			# If the line starts with a ">" sign, it is a read identifier line
			if lines.startswith(">"):
				
				# If the end of the line is not a "]", we need to concatenate this line with the next line
				if lines[len(lines) - 1] != "]":
					
					# Add the current line to the empty string s_cont and set the boolean variable to True
					s_cont = s_cont + lines
					cont = True
					
					# Immediately move forward to the next iteration to get the next line
					continue
				
				# If the end of line is a "]", the species name can be extracted through this line
				else:
					# Get the species name given the read identifier line
					name = get_species(lines)
					
					# Since the species name is present as a significant hit with a key enzyme in CpG methylation system,
					# its name is removed from the original list and added to the list of species that might have a 
					# methylation system
					add_or_remove_species(name, species_with_ms, species)
					
		# If there is a need to concatenate two lines
		else:
			
			# Append a whitespace followed by the current line to s_cont, which contains the previous line content
			s_cont = s_cont + " " + lines
			
			# Get the species name from the concatenated lines
			spec = get_species(s_cont)
			
			# Since the species name is present as a significant hit with a key enzyme in CpG methylation system,
			# its name is removed from the original list and added to the list of species that might have a 
			# methylation system
			add_or_remove_species(spec, species_with_ms, species)
			
			# Reset both s_cont and cont (s_cont--> empty string | cont --> False) so that they can be used again. 
			s_cont = ""
			cont = False

# Prints out results 
print("Here is a list of hexapoda species that might have a methylation system: " + "\n")
print(species_with_ms)
print()
print("Here is a list of hexapoda species that has no evidence of a methylation system: " + "\n")
print(species)
print()
print("Out of " + str(total_species) + " hexapoda species, " + str(len(species)) + " of them seem not to have a methylation system")
		
