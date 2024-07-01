#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script concatenates sequences from a set of fasta files (input as a list) into a new concatenation fasta.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)

import os
import re
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"Bio" : "https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython"}
for nstlobject, link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

from Bio.SeqIO.FastaIO import SimpleFastaParser

print('#Script: concatenation.py')
print('#Version: v20240701')
print('#Usage: python concatenation.py <dataset_names> <output_file>')
print("#<dataset_names> must be a list (one per line) of the files/alignments that will be concatenated. (required)")
print('#<output_file> must be the name of the output concatenation FASTA file (can contain alphanumeric characters, underscores, and dots). (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

#Create empty list for unique header accessions (those should be assembly accessions)
headers_accessions = []

#Checkpoint for number of arguments
if len(sys.argv) == 3:
    print ('Two arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

#Checkpoint for the existence of the dataset_names file
if os.path.isfile(sys.argv[1]) == True:
    print('Dataset names file found. Proceeding.')
else:
    print('Dataset names file not found. Exiting.')
    sys.exit(1)

#Checkpoint for disallowed characters in output file name.
#TODO: Maybe remove this character requirement in the output filename? It's not adding anything nor will it prevent doggo_zoomies or iqtree from running. Otherwise remove dots as an allowed character and automatically use .concatenation as the extension.
if re.match(r'^[A-Za-z0-9_\.]+$', sys.argv[2]):
    print ('Concatenation name is valid. Proceeding.')
else:
    print ('Concatenation name is invalid. Exiting.')
    sys.exit(1)

#Checkpoint for the existence of the files in the dataset list.
with open(sys.argv[1]) as f:
	datasets = [line.rstrip() for line in f]
for filename in datasets:
    if not os.path.isfile(filename):
        print('Dataset ' + filename + ' not found. Exiting.')
        sys.exit(1)
print('All dataset files found. Proceeding.')

#Remove any previous output files with the same name.
print ('Removing files with names identical to the output.')
removal = ('rm -r ' + sys.argv[2] + ' 2> /dev/null')
os.system(removal)

#Open list with filenames to be concatenated.
#TODO: Check if files in list are correctly formatted fasta?
#filenames = [sys.argv[1]
print('Extracting sequence accessions.')
filedata = {filename: open(filename, 'r') for filename in datasets}
#Second argument is concatenation output file.
concat_out = open(sys.argv[2], 'a')
#Create a list with all the different accessions.
for file in filedata.values():
    #For each sequence in each dataset, if the assembly is not in list of unique accessions, add it.
	for title1, seq1 in SimpleFastaParser(file):
		x = title1.split(' ')
		if x[0] not in headers_accessions:
			headers_accessions.append(x[0])

print('Creating concatenation.')
#Loop through the gene files and concatanate with the accession list.
#For each unique accession, create an empty concatenated sequence.
for item in headers_accessions:
	concat_seq = ''
#For each of the datasets, create a null variable D that determines presence or absence of sequence.
	filedata = {filename: open(filename, 'r') for filename in datasets}
	for file in filedata.values():
		D = 0
#For each of the sequences, determine their length (N), extract accessions.
		for title, seq in SimpleFastaParser(file):
			N = len(seq)
			y = title.split(' ')
#If you find a sequence matching the unique accession currently checked, add the sequence to the concatenated sequence. D increase by 1.
			if item == y[0]:
				concat_seq = concat_seq + seq
				D+=1
#If no sequence is found, D is 0, add to the sequence dashes (indels) equal to the length.
		if D == 0:
			concat_seq = concat_seq + '-'*N
#Write to the output file the concatenated sequence with the accession as header.
	concat_out.write('>' + item + '\n' + concat_seq + '\n')
concat_out.close()

print('All done!')
