#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script does a dry run of pickgenomes, reporting from which taxa and how many genomes will be picked, and the total number.

#Dependencies
#NONE

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

from collections import Counter
import os
import sys

print('#Script: pickgenomes_dry.py')
print('#Version: v20240713')
print('#Usage: python pickgenomes_dry.py <input_tsv> <tax_level> <tax_resolution> <number> <ignore_list>')
print('#<input_tsv> must be a tab-delimited file as per GTDB\'s metadata files for Bacteria or Archaea. (required)')
print('#<tax_level> must be the highest taxonomic level for genomes to be selected as presented in GTDB taxonomy strings e.g. p__Asgardarchaeota (i.e., get genomes from within <tax_level>). (required)')
print('#<tax_resolution> must be a taxonomic level lower than <tax_level>. For each <tax_resolution> in <tax_level>, the script picks <number> genomes (or as many as available). It must be p (phylum), c (class), o (order), f (family), g (genus), or all. "all" picks all genomes for <tax_level>. (required)')
print('#<number> must be the number of genomes to be picked per <tax_resolution> or "all". If <tax_resolution> is "all", <number> must also be "all". (required)')
print('#<ignore_list> must be a text file containing genome assembly accessions (once per line, versionless, e.g. GCA_011362025) that should be not included in the pickgenomes process (optional)')
print('#For more information refer to the comments in the script and/or the Github page.')

resolutions = list()
ignore_list = list()
all_res_count = 0

# checkpoint for number of arguments
if len(sys.argv) == 5 or len(sys.argv) == 6:
    print(str((len(sys.argv)-1)) + ' arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

file_stem = (str(sys.argv[2]) + '_' + str(sys.argv[3]) + '_' + str(sys.argv[4]))

# checkpoint for input file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Input file found. Proceeding.')
else:
    print('Input file not found. Exiting.' + '\n')
    sys.exit(1)

#Check if input file is parsed GTDB metadata.
with open(sys.argv[1], 'r') as parsedcheck:
    for line in parsedcheck:
        x = line.split('\t')
        if x[18]=='t':
            pass
        else:
            print('Input file does not contain parsed GTDB metadata. Exiting.')
            sys.exit(1)
print('Input file contains parsed GTDB metadata. Proceeding.')

# checkpoint for taxonomic level
with open(sys.argv[1], 'r') as file:
    content = file.read()
    if str(sys.argv[2]+';') in content:
        print('Taxonomic level is valid. Proceeding.')
    else:
        print('Taxonomic level is invalid. Check spelling. Exiting.' + '\n')
        sys.exit(1)

# checkpoint for taxonomic resolution
tax_res = ['p', 'c', 'o', 'f', 'g']
if sys.argv[3] in tax_res:
    print('Taxonomic resolution is valid. Proceeding.')
elif sys.argv[3] == 'all':
    print('Taxonomic resolution is valid (all). Proceeding.')
else:
    print('Taxonomic resolution is invalid or missing. Check if it is one of: "p", "c", "o", "f", "g", or "all". Exiting.' + '\n')
    sys.exit(1)

# checkpoint for resolution and number of genomes both being "all"
if sys.argv[3] == 'all' and sys.argv[4] != 'all':
    print('When <tax_resolution> is "all", <number> must also be "all". Exiting.' + '\n')
    sys.exit(1)
elif sys.argv[3] != 'all' and sys.argv[4] == 'all':
    print('When <number> is "all", <tax_resolution> must also be "all". Exiting.' + '\n')
    sys.exit(1)

# checkpoint for number of genomes to be picked
flag = True
try:
    int(sys.argv[4])
except:
    flag = False
if flag and int(sys.argv[4])>0:
    print('Number given is valid. Proceeding.')
elif sys.argv[4] == 'all':
    print('Number given is valid (all). Proceeding.')
else:
    print('Number given is invalid or missing. Number must be a natural number or both <tax_resolution> and <number> must be "all". Exiting.' + '\n')
    sys.exit(1)

if len(sys.argv) == 6: ## this is to make the ignore_list argument optional
    check_ignore = os.path.isfile(sys.argv[5])
    if check_ignore == True:
        print('Ignore list file found. Proceeding.')
        with open(sys.argv[5], 'r') as ignore_text:
            for line in ignore_text:
                x=line.strip()
                ignore_list.append(x)
    else:
        print('Ignore list file not found. Exiting.')
        sys.exit(1)
else:
        print('No ignore list specified. Proceeding.')

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r ' + file_stem + '.dry 2> /dev/null')
os.system(removal)

print('Picking genomes.')

with open(sys.argv[1], 'r') as taxa_records:
    for line in taxa_records:
        x = line.split('\t')
        if x[0] != 'accession': #Ignore the headers line
            y = x[19].split(';')
            x[0] = x[0][:5] + 'A' + x[0][6:] if x[0][5] == 'F' else x[0]
            if x[0][3:-2] not in ignore_list and sys.argv[2] in y: ## this is the step where we exclude all the assemblies from the ignore_list
                if sys.argv[3] == str(y[1])[0]: #if <tax_resolutions> equals the first letter of y element 1 (aka equals "p" for phylum) add it to the resolutions list (will have duplicates).
                    resolutions.append(y[1])
                elif sys.argv[3] == str(y[2])[0]: #ditto for class
                    resolutions.append(y[2])
                elif sys.argv[3] == str(y[3])[0]: #ditto for order
                    resolutions.append(y[3])
                elif sys.argv[3] == str(y[4])[0]: #ditto for family
                    resolutions.append(y[4])
                elif sys.argv[3] == str(y[5])[0]: #ditto for genus
                    resolutions.append(y[5])
                elif sys.argv[3] == 'all':
                    all_res_count += 1

resolution_counts = Counter(resolutions)
total_genomes_tbd = 0 # total number genomes to be downloaded

with open(file_stem + '.dry', 'w') as outdry:
    if sys.argv[4] == 'all':
        outdry.write("Number of " + sys.argv[2] + " genomes to be downloaded: " + str(all_res_count))
    else:
        for resolution, count in resolution_counts.items(): # resolutions has multiplicates
            if count <= int(sys.argv[4]):
                total_genomes_tbd = total_genomes_tbd + count
                outdry.write(f"Number of {resolution} genomes to be downloaded: {count}" + "\n")
            else:
                total_genomes_tbd = total_genomes_tbd + int(sys.argv[4])
                outdry.write(f"Number of {resolution} genomes to be downloaded: " + sys.argv[4] + "\n")
        outdry.write("Total number of genomes to be downloaded: " + str(total_genomes_tbd))

print('All done!')
