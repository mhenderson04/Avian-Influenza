import os
from Bio import SeqIO
import glob

def main():

    ### Parse the file for the gene of interest ###
    while True:
        # Call file to parse
        fileInput = input("Enter file name. If you are done using the program, enter 'exit': ")
        file = fileInput
        if file == 'exit':
            break

        # Enter gene to search for
        geneInput = input("Enter gene name: ")
        gene = geneInput

        # Pass varibales to parse function
        parse_fasta_by_description(file, gene)

def parse_fasta_by_description(input_file, description_pattern):
    # Open the input FASTA file
    with open(input_file, "r") as infile:
        records = SeqIO.parse(infile, "fasta")
        
        # Initialize a counter for file numbering
        file_counter = 1
        current_output_file = None
        current_output_handle = None
        
        # Iterate through each record
        for record in records:
            # Check if the description contains the pattern
            if description_pattern in record.description:
                # Close the previous file if one was opened
                if current_output_handle:
                    current_output_handle.close()
                
                # Create a new output file for this record if necessary
                output_file = f"{description_pattern}.fasta"
                current_output_handle = open(output_file, "a")
                
                # Write the record to the current output file
                SeqIO.write(record, current_output_handle, "fasta")
                
                # Increment the file counter if you wish to write a new file per record
                file_counter += 1
        
        # Close the final output file if opened
        if current_output_handle:
            current_output_handle.close()

main()

