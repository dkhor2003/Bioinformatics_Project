#!/opt/rit/spack-app/linux-rhel7-x86_64/gcc-10.2.0/python-3.8.8-ucekvffll3knventjp2555etn5ee7jqu/bin/python

# Import Bio.Entrez to download protein sequences
from Bio import Entrez
import sys

# Specify the output file to download the protein sequences in FASTA format
outfile = sys.argv[1] 

# User email
Entrez.email = 'dkhor@iastate.edu'

# Searches for the accession IDs for the protein sequences desired in the NCBI protein database 
entrez_handle = Entrez.esearch(db = 'protein', term = '("Apis mellifera"[Organism] AND DNMT1a[Gene Name]) OR ("Apis mellifera"[Organism] AND DNMT3[Gene Name]) AND refseq[filter]', idtype = "acc")
entry = Entrez.read(entrez_handle)

# Get the list of accession IDs from the parsed entrez_handle
ids = list(entry["IdList"])

# Retrieves the protein sequence records by using the accession ID list
entrez_handle2 = Entrez.efetch(db = 'protein', id = ids, rettype = 'fasta', retmode = 'txt')
entry2 = entrez_handle2.read()

# Closes both handles used
entrez_handle.close()
entrez_handle2.close()


# Write the protein sequences into the specified outfile
with open(outfile, "w") as file:
	file.write(entry2)

# Open the outfile that has been writen and print out the sequence names (FASTA identifiers)
with open(outfile, "r") as file2:
	for line in file2:
		if line.startswith(">"):
			print(line.lstrip(">"))
