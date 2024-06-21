#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script recodes sequences in amino acid fasta files to different reduced amino acid alphabets.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.
#NOTE 2: We chose to do all recoding to numbers, which works for IQ-TREE, even the 4-state alphabets.

#Dependencies
#NONE

import os
import sys

print('#Script: aarecode.py')
print('#Version: v20240621')
print('#Usage: python aarecode.py <input_faa> <output_faa> <alphabet>')
print('#<input_faa> must be the input protein FASTA file using the complete 20-state amino-acid alphabet. Unknown states are assumed to be X. (required)')
print('#<output_faa> must be the name of the output recoded FASTA file. (required)')
print('#<alphabet> must be the reduced amino acid alphabet. It must be D6 (Dayhoff 6-state) or SR4 (Susko-Roger 4-state). (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 4:
    print ('Three arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

# checkpoint for input file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Input file found. Proceeding.')
else:
    print('Input file not found. Exiting.')
    sys.exit(1)
#TODO: Add a check that the file is really in fasta format?

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r ' + sys.argv[2] + ' 2> /dev/null')
os.system(removal)

# Checkpoint for the reduced aa alphabet
schemes = ('D6', 'SR4')
if sys.argv[3] in schemes:
    print('Reduced amino acid alphabet is valid. Proceeding.')
else:
    print('Reduced amino acid alphabet is not valid. Exiting.')
    sys.exit(1)

# Define the reduced amino acid alphabets here
def recode_sequence(sequence, red_aa):
    if red_aa == "D6":
        translation_dict = {
            "A": "0", "G": "0", "P": "0", "S": "0", "T": "0",
            "D": "1", "E": "1", "N": "1", "Q": "1",
            "H": "2", "K": "2", "R": "2",
            "F": "3", "Y": "3", "W": "3",
            "I": "4", "L": "4", "M": "4", "V": "4",
            "C": "5"
        }
        new_sequence = ""
        for aa in sequence:
            if aa == 'X':
                new_sequence += '?'
            else:
                new_sequence += translation_dict.get(aa, aa)
        return new_sequence
    elif red_aa == "SR4":
        translation_dict = {
            "A": "0", "G": "0", "N": "0", "P": "0", "S": "0", "T": "0",
            "F": "1", "I": "1", "L": "1", "M": "1", "V": "1",
            "C": "2", "H": "2", "W": "2", "Y": "2",
            "D": "3", "E": "3", "K": "3", "Q": "3", "R": "3"
        }
        new_sequence = ""
        for aa in sequence:
            if aa == 'X':
                new_sequence += '?'
            else:
                new_sequence += translation_dict.get(aa, aa)
        return new_sequence
    # Add more alphabets here as needed
    else:
        return sequence

print('Recoding sequences.')
with open(sys.argv[1], "r") as goat, open(sys.argv[2], "w") as sheep:
    for line in goat:
        if line.startswith(">"):
            sheep.write(line)  # Keep header lines
        else:
            recoded_sequence = recode_sequence(line.strip(), sys.argv[3])
            sheep.write(recoded_sequence + "\n")  # Perform recoding on sequence lines

print('All done!')
