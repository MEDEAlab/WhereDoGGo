#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script removes all sequences with the same accession number from a fasta file, outputting the rest in a new fasta, and noting the accessions removed and number of multiplicates in a log.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.
#NOTE 2: This is the WhereDoGGo? version of the script that works on all files in the working directory with a given extension.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)

import os
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"Bio" : "https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

from Bio import SeqIO

print('#Script: removemultiples.py')
print('#Version: v20241212')
print('#Usage: python removemultiples.py <datasets> <input_ext> <output_ext> <output_log>')
print('#<datasets> must be the directory containing the FASTA files with <input_ext>. (trailing slash optional) (required)')
print("#<input_ext> must be the extension of the FASTA files in <datasets> that will be checked for sequences with the same accession. The stem of each file is retained for the output and log files. (leading dot optional) (required)")
print('#<output_ext> must be the extension of the created FASTA files where the sequences without multiples will be written. (leading dot optional) (required)')
print('#<output_log> must be the name of the tab-delimited log file that will contain the accessions with multiples and the number of times each was found. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

#Checkpoint for number of arguments
if len(sys.argv) == 5:
    print ('Four arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

#Check if the extensions start with a dot, otherwise add them.
input_ext = sys.argv[2]
output_ext = sys.argv[3]
if input_ext.startswith('.') == False: #This is not absolutely necessary, since os.path.splitext will detect the extension anyway. It's more of a precaution against double dots.
    input_ext = str('.' + input_ext)
if output_ext.startswith('.') == False:
    output_ext = str('.' + output_ext)

#Checkpoint for datasets directory existence and trailing slash. Convert to abspath to make sure there are no issues when called through doggo_sniff.
if os.path.exists(sys.argv[1]) == True:
    print ('Datasets directory found. Proceeding.')
    datasetsdir = os.path.abspath(sys.argv[1])
    datasetsdir = os.path.join(datasetsdir, '')
else:
    print ('Datasets directory not found. Exiting.')
    sys.exit(1)

#Check if files with a given extension exist in the datasets directory and create a list of them.
filenames = []
for fname in os.listdir(datasetsdir):
    if fname.endswith(input_ext):
        fname = os.path.join(datasetsdir, fname)
        filenames.append(fname)
if len(filenames) > 0:
    print('File(s) with the input extension found in the datasets directory. Proceeding.')
else:
    print('No files with the input extension found in the datasets directory. Exiting.')
    sys.exit(1)

#Remove any previous output files with the same name.
print ('Removing files with names identical to the output.')
removal = ('rm -r *' + output_ext + ' ' + sys.argv[4] + ' 2> /dev/null')
os.system(removal)

print ('Writing log and output files.')
output_log = open(sys.argv[4], "w")
for fname2 in filenames:
    #Separate the fasta file stem to use for the output and log files.
    filestem = str(os.path.basename(fname2).split(os.extsep, 1)[0])
    filestemplusinputext = str(os.path.basename(fname2))
    #Create an empty accessions list and an empty removals list. Populate the accessions list with all accessions in the current fasta file.
    accessions = []
    removals = []
    output_log.write(filestemplusinputext + '\n')
    for acc in SeqIO.parse(fname2, 'fasta'):
        accessions.append(acc.id)
    #For each accession, if it appears more than once, add it to the removals set. This will also remove the multiplicates from the list.
    for accrem in accessions:
        howmany = accessions.count(accrem)
        if howmany > 1 and accrem not in removals:
            removals.append(accrem)
            howmany = str(howmany) #NOTE: Since python can't write an integer like this, the number get assigned to a variable and the variable is converted to a string.
            output_log.write(accrem + '\t' + howmany + '\n')
    output_log.write('//' + '\n')
    #Open the output fasta file. Write inside any sequences that are NOT in the set with multiplicates. Then close it.
    #TODO: Parsing the fasta file twice is inefficient. Need to find a workaround.
    leftover = (acckeep for acckeep in SeqIO.parse(fname2, "fasta") if acckeep.id not in removals)
    output_seq = open(str(filestem + output_ext), "w")
    SeqIO.write(leftover, output_seq, "fasta")
    output_seq.close()
output_log.close()

print('All done!')
