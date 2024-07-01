#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script outputs a tab-delimited file with the taxonomic distribution of sequences in each FASTA file with a given extension (assembly accession, name, tab-delimited sequence accessions), based on a tab-delimited assembliesnames file.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.
#NOTE 2: This is the WhereDoGGo? version of the script that works on all files in the working directory with a given extension.

#Dependencies
#1) pandas (https://github.com/pandas-dev/pandas)

import csv
import os
import re
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"pandas" : "https://github.com/pandas-dev/pandas"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

import pandas as pd

print('#Script: fasta2distribution.py')
print('#Version: v20240701')
print('#Usage: python fasta2distribution.py <input_ext> <distribution_file> <assembliesnames_file>')
print("#<input_ext> must be the filename extension of the input FASTA files. It is assumed that the header will start with the protein accessions, followed by the assembly accession, and then anything else, separated by a space. (required)")
print('#<distribution_file> must be the name of the output file that will contain the taxonomic distribution of the sequences in the input FASTA files. (required)')
print('#<assembliesnames_file> must be a tab-delimited file with all assembly-name pairs. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 4:
    print ('Three arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

input_ext = sys.argv[1]
if input_ext.startswith('.') == False: #This is not absolutely necessary, since os.path.splitext will detect the extension anyway. It's more of a precaution against double dots.
    input_ext = str('.' + input_ext)

#Check if files with a given extension exist in the working directory and create a list of them.
filenames = []
for fname in os.listdir('.'):
    if fname.endswith(input_ext):
        filenames.append(fname)
if len(filenames) > 0:
    print('File(s) with the input extension found in the working directory. Proceeding.')
else:
    print('No files with the input extension found in the working directory. Exiting.')
    sys.exit(1)

if os.path.isfile(sys.argv[3]) == True:
    print('Assembliesnames file found. Proceeding.')
else:
    print('Assembliesnames file not found. Exiting.')
    sys.exit(1)

#Remove any previous output files with the same name.
print ('Removing files with names identical to the output.')
removal = ('rm -r *.distro ' + sys.argv[2] + ' 2> /dev/null')
os.system(removal)

print ('Creating individual distribution files.')
assemblies_names_dict = {}  # dictionary { assembly : name }
with open(sys.argv[3], 'r') as file:
    for line in file:
        assembly, name = line.strip().split('\t')
        assemblies_names_dict[assembly] = name

for fname2 in filenames:
    #separate the fasta file stem to use for the output and log files.
    filestem = os.path.splitext(fname2)[0]
    with open(fname2, 'r') as infile, open(str(filestem + '.distro'), 'w') as outfile:
        for assembly, name in assemblies_names_dict.items():
            distribution = str(assembly + '\t' + name)
            assembly_found = False  # Flag to check if assembly is found in any line
            # this resets the file cursor to the beginning of the file for each assembly
            infile.seek(0)
            for line in infile:
                if assembly in line:
                    prot_acc_match = re.search(r'>(\S+)', line)  # protein accession regex
                    if prot_acc_match:
                        prot_acc = prot_acc_match.group(1)
                        distribution = distribution + '\t' + prot_acc
                    assembly_found = True  # Set the flag to True if assembly is found
            outfile.write(distribution + '\n')
        if not assembly_found:
            # if assembly is not found in any line, write the default distribution
            outfile.write(distribution + '\n')

# add headers to the distribution files
for file_path in os.listdir('.'):
    # Read the file and split rows
    if file_path.endswith('.distro'):
        with open(file_path, 'r') as file:
            rows = [line.strip().split('\t') for line in file]

        # Find the row with the maximum number of substrings, needed to add NAs afterwards
        max_substrings_row = max(rows, key=len)

        filestem2 = os.path.splitext(os.path.basename(file_path))[0]

        headers = ['taxon_name', 'assembly'] + [filestem2] * (len(max_substrings_row) - 2)

        # Write the headers and processed rows back to the file
        with open(file_path, 'w') as file:
            file.write('\t'.join(headers) + '\n')
            for row in rows:
                file.write('\t'.join(row) + '\n')

# add NA to empty cells
for file_path in os.listdir('.'):
    # Read the file and split rows
    if file_path.endswith('.distro'):
        df = pd.read_csv(file_path, sep='\t', header=0, index_col=False, na_values=[''], dtype=str)
        df.fillna('NA', inplace=True)
        df.to_csv(file_path, sep='\t', index=False)

print ('Creating combined distribution file.')
# make concatenated distribution file
concatenated_matrix = pd.DataFrame()
files = [file for file in os.listdir('.') if file.endswith('.distro')]

for file in files:
    df = pd.read_csv(file, sep='\t', dtype=str)

    # Extract the first two columns (assuming they are named 'taxon_name' and 'assembly')
    first_two_columns = df[['taxon_name', 'assembly']]

    # Extract all columns after the second one (headers)
    remaining_columns = df.iloc[:, 2:]

    # Concatenate the DataFrames along columns
    concatenated_matrix = pd.concat([concatenated_matrix, remaining_columns], axis=1)

# Concatenate the first two columns to the concatenated matrix
final_matrix = pd.concat([first_two_columns, concatenated_matrix], axis=1)

# Remove '.XX' from the headers
final_matrix.columns = [col.split('.')[0] for col in final_matrix.columns]

# Write the final matrix to a new tab-separated file
final_matrix.to_csv(sys.argv[2], sep='\t', index=False)

print ('Removing individual distribution files.')
# Remove individual distribution files
for file_path in files:
    os.remove(file_path)

print('All done!')
