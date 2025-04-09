import os

def main():
    dirName = input("Enter directory to store file: ")
    outDir = dirName
    fileName = input("Name the output file: ")
    outFile = fileName
    merge_fasta_files(outDir, outFile)


def merge_fasta_files(output_directory, final_output_file):
    with open(final_output_file, 'w') as outfile:
        for file_name in os.listdir(output_directory):
            if file_name.endswith('.fasta'):
                file_path = os.path.join(output_directory, file_name)
                
                with open(file_path, 'r') as infile:
                    outfile.write(infile.read() + "\n\n")  # Add newline between files
                
                print(f"Merged: {file_path}")
main()
