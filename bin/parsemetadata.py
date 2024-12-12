#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script will parse a file containing GTDB metadata and output only the lines containing the representative genome for each species cluster.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.

#Dependencies
#NONE

import os
import sys

print('#Script: parsemetadata.py')
print('#Version: v20241212')
print('#Usage: python parsemetadata.py <input_file> <output_file>')
print('#<input_file> must be tab-delimited GTDB metadata. (required)')
print('#<output_file> must be the name of the output file that will contain the representative genome of each species cluster. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 3:
    print ('Two arguments found. Proceeding.')
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
#TODO: Check if the input file is what's expected, at least if it's tab-delimited.

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r ' + sys.argv[2] + ' 2> /dev/null')
os.system(removal)

#Keep lines marked with a "t" in the 16th field, denoting representative genome.
print('Parsing GTDB metadata')
with open(sys.argv[1], 'r') as metadata, open(sys.argv[2], 'w') as parsed:
    for line in metadata:
        x = line.split('\t')
        if x[18]=='t':
            parsed.write(line)

print('All done!')
