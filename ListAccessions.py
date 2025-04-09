from Bio import SeqIO

# Prompt for file name(s)
fileInput = input("Enter file name: ")
fastaFile = fileInput

# Parse the file and list accessions
with open(fastaFile, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        accession = record.id.split()[0]
        print(accession)
