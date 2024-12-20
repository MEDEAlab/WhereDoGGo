#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script checks amino acid sequences in a fasta file, converting lowercase to uppercase and converting non-standard amino acids and other symbols to "X" for unknown state, to avoid weird behavior by other programs.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.

#Dependencies
#NONE

import os
import re
import sys

print('#Script: fixdualsupports.py')
print('#Version: v20241212')
print('#Usage: python fixdualsupports.py <input_tree> <output_tree>')
print('#<input_tree> must be the input phylogenetic tree in Newick format, with two slash-separated support values of which the first must be ultrafast bootstraps and the second aLRT-SH. (required)')
print('#<output_tree> must be the name of the output Newick tree file. Strongly supported branches (UFBOOT>=95 & aLRT-SH>=0.80) are assigned a value of 1 and weakly supported branches are assigned a value of 0. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

#Check if the correct number of arguments is given
if len(sys.argv) == 3:
    print('Two arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

#checkpoint for input tree file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Input tree file found. Proceeding.')
else:
    print('Input tree file not found. Exiting.')
    sys.exit(1)
#TODO: Check for Newick format?

#Remove any previous output files with the same name.
print('Removing files with names identical to the output.')
removal = ('rm -r ' + sys.argv[2] + ' 2> /dev/null')
os.system(removal)

print('Fixing branch supports.')
#Argument 2: Define the output tree file
outtree = open(sys.argv[2], 'a')

#Argument 1: Read the input tree file
with open(sys.argv[1], 'r') as intree:
    #Read all lines from the input tree
    lines = intree.readlines()
    #Read the entire input tree into a variable
    filedata = "".join(lines)

    #Initialize a list to store supports found in the lines
    supports = []

    #Use regular expressions to find supports
    current_supports = re.findall(r"(\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)(?=:)", filedata)
    supports.extend(current_supports)

    #Process each support found
    for support in supports:
        #Extract numerator and denominator from the tuple
        numerator, denominator = map(float, support)
        #Check if the current file data contains the support value
        support_str = f"{support[0]}/{support[1]}"
        if support_str in filedata:
            if numerator >= 95 and denominator >= 80:
                #Replace supports of ultra-fast bootstrap that are >=95 and aLRT-SH >= 80 with 1
                filedata = filedata.replace(support_str, "1", 1)
            else:
                #Replace supports of ultra-fast bootstrap that are <95 or aLRT-SH < 80 with 0
                filedata = filedata.replace(support_str, "0", 1)

    #Write the modified data to the output tree file
    outtree.write(filedata)

#Close the output tree file
outtree.close()

print('All done!')
