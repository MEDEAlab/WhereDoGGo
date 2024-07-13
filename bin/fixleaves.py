#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script changes the leaf names that correspond to assembly accessions in a Newick tree file, to the corresponding lines in an assembliesnames file (tah-delimited becomes space-delimited).

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#1) ETE3 (https://anaconda.org/etetoolkit/ete3 or https://pypi.org/project/ete3/)

import os
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"ete3": "https://anaconda.org/etetoolkit/ete3 or https://pypi.org/project/ete3/"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

from ete3 import Tree

print('#Script: fixleaves.py')
print('#Version: v20240713')
print('#Usage: python fixleaves.py <input_tree> <assembliesnames_file> <output_tree>')
print('#<input_tree> must be the input phylogenetic tree in Newick format. (required)')
print('#<input_names> must be a file containing the names of the leaves that will be used for renaming. The current leaf name must be contained in the new names. Tabs as in .assembliesnames are converted to spaces. (required)')
print('#<output_tree> must be the name of the output Newick tree file. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 4:
    print ('Three arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

# checkpoint for input tree file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Input tree file found. Proceeding.')
else:
    print('Input tree file not found. Exiting.')
    sys.exit(1)
#TODO: Check for Newick format?

# checkpoint for assembliesnames file
check_file = os.path.isfile(sys.argv[2])
if check_file == True:
    print('Assembliesnames file found. Proceeding.')
else:
    print('Assembliesnames file not found. Exiting.')
    sys.exit(1)
#TODO: Check for formatting?

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r ' + sys.argv[3] + ' 2> /dev/null')
os.system(removal)

print('Converting leaf names.')
# Open the output tree file for writing
with open(sys.argv[3], "w") as output_tree_file:
    # Open the assembliesnames file and read its content
    with open(sys.argv[2], "r") as assembliesnames_file:
        species_assembly_dict = {}
        for line in assembliesnames_file:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                assembly_accession, species_name = parts
                species_assembly_dict[assembly_accession] = species_name
        # Open the input tree file and process each Newick tree
        with open(sys.argv[1], "r") as input_tree_file:
            newick_trees = input_tree_file.read().split(';')
            for newick_tree in newick_trees:
                if not newick_tree.strip():
                    continue  # Skip empty entries or entries with just whitespace
                try:
                    # Parse the Newick tree
                    tree = Tree(newick_tree + ";")
                    # Update the leaf names in the tree
                    for leaf in tree.get_leaves():
                        accession = leaf.name
                        species_name = species_assembly_dict.get(accession, "")
                        if species_name:
                            leaf.name = f"{accession} {species_name}"
                    # Write the updated tree to the output file
                    output_tree_file.write(tree.write(format=0) + "\n")
                except:
                    print(f"Error processing the following tree: {newick_tree}")
                    sys.exit(1)

print('All done!')
