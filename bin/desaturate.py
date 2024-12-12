#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script removes positions fast-evolving from a fasta alignment using a file with site-specific rates, to produce a number of desaturated alignments based on a given number of quantiles.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)
#2) pandas (https://github.com/pandas-dev/pandas)
#3) NumPy (https://numpy.org/install/)

import math
import os
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"Bio":"https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython",
                        "numpy":"https://numpy.org/install/",
                        "pandas":"https://github.com/pandas-dev/pandas"}
for nstlobject,links in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

from Bio import AlignIO
import numpy as np
import pandas as pd

print('#Script: desaturate.py')
print('#Version: v20241212')
print('#Usage: python desaturate.py <rates> <quantiles> <alignment>')
print('#<input_rates> must be a file containing site-specific rates for <input_alignment> produced by IQ-TREE\'s --rate or --mlrate options. (required)')
print('#<input_quantiles> must be the number (integer) of quantiles of retained positions. The original alignment, the empty alignment, and anything with more than 90% of positions removed are omitted. For example, for <input_quantiles> = 20, 18 alignments will be produced from 5% of positions retained, up to 90%. (required)')
print('#<input_alignment> must be the input alignment in FASTA format. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Create an empty dictionary to store site rates
site_rate = dict()
# Create an empty list to store rates
rates = list()

# Check if the correct number of arguments is given
if len(sys.argv) == 4:
    print ('Three arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

# checkpoint for input tree file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Rates file found. Proceeding.')
else:
    print('Rates file not found. Exiting.')
    sys.exit(1)
#TODO: Check for formatting?

# Check if the number of desired quantiles is a valid integer between 0 and 100
if (isinstance(int(sys.argv[2]), int) and 0 < int(sys.argv[2]) < 100):
    print("Number of quantiles is valid. Proceeding.")
else:
    print('Number of quantiles is not valid. Exiting.')
    sys.exit(1)

# Check if the input alignment file exists
check_file = os.path.isfile(sys.argv[3])
if check_file == True:
    print('Alignment file found. Proceeding.')
else:
    print('Alignment file not found. Exiting.')
    sys.exit(1)
#TODO: Check if alignment is in fasta format?

#Isolate the stem of the input alignment file to be used as the stem for the output alignments. This means that the script will always output in the working directory.
outfile_stem = str(os.path.basename(sys.argv[3]).split(os.extsep, 1)[0])

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r ' + outfile_stem + '_removed*.faa 2> /dev/null')
os.system(removal)

# Read rates from the input rates file and create a dictionary
with open(sys.argv[1], 'r') as rate_file:
    for line in rate_file:
        if not line.startswith('#') and not line.startswith('Site'):
            x = line.split('\t')
            site_rate.update({x[0]: x[1].strip()})  # Create a dictionary with key = site, value = rate

# Extract rates from the dictionary and convert them to floats
for key, value in site_rate.items():
    rates.append(float(value))

# Create a Pandas Series from the rates
ser = pd.Series(rates)

# Perform quantile-based discretization and get the bin edges
results, bins_edges = pd.qcut(x=ser, q=int(sys.argv[2]), retbins=True, duplicates='drop')

# Initialize lists for sites to be deleted and bin edges
sites_to_be_deleted = list()
bins_edges_list = list()

# Convert bin edges to a list
for bin_edge in bins_edges:
    bins_edges_list.append(float(bin_edge))

# Sort the bin edges in descending order
bins_edges_list_reversed = sorted(bins_edges_list, reverse=True)

# Calculate the increment for site removal and set the initial removed count
increment = 100 / int(sys.argv[2])
increment = math.ceil(increment) # Round up if necessary, no decimals
removed_count = 0 - int(increment)  # Negative value to start desaturation from 0

print('Creating desaturated datasets.')
# Loop through bin edges for desaturation
for bin_edge in bins_edges_list_reversed:
    removed_count = removed_count + int(increment)
    if removed_count != 0 and removed_count <= 90:  # Don't remove 0 sites and don't remove more than 90% of the sites
        # Create an output file with the removed count in the name
        outfile = open((outfile_stem + '_removed' + str(removed_count)) + '.faa', 'a')
        # Read the input alignment in FASTA format
        alignment = AlignIO.read(sys.argv[3], "fasta")

        # Identify sites to be deleted based on the rate
        for key, value in site_rate.items():
            if float(value) > float(bin_edge):
                new_key = int(key) - 1  # Subtract 1 from the site number because the numpy array starts from 0
                sites_to_be_deleted.append(int(new_key))

        # Sort the sites to be deleted in descending order
        sites_to_be_deleted_reversed = sorted(sites_to_be_deleted, reverse=True)

        # Remove the selected sites from the alignment
        for stbd in sites_to_be_deleted_reversed:
            alignment = alignment[:, :stbd] + alignment[:, (int(stbd) + 1):]

        # Write the modified alignment to the output file
        AlignIO.write(alignment, outfile, "fasta")

        # Clear the list of sites to be deleted for the next iteration
        sites_to_be_deleted.clear()

# Close the output file.
outfile.close()

print('All done!')
