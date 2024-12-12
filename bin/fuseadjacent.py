#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2023 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script fuses fragmented ORFs (determined from successive accession numbers) from a fasta file into a new fasta file.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.
#NOTE 2: The script assumes that the FASTA headers were produced by Pyrodigal, as used in doggo_sniff.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)

import os
import re
import statistics
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
from Bio.Seq import Seq

print('#Script: fuseadjacent.py')
print('#Version: v20241212')
print('#Usage: python fuseadjacent.py <datasets> <input_ext> <output_ext> <output_log>')
print('#<datasets> must be the directory containing the FASTA files with <input_ext>. (trailing slash optional) (required)')
print('#<input_ext> must be the extension of the alignment FASTA files in <datasets> that will be checked for adjacent fragmented sequences. Temporary files with unfused and fused sequences in single-line FASTA (aka nice FASTA) will be created. The stem of each file is retained for the temporary and output files. (leading dot optional). (required)')
print('#<output_ext> must be the extension of the FASTA files where the fused and unfused sequences sequences will be written in single-line format. (leading dot optional) (required)')
print('#<output_log> must be the name of the output log file that will contain the fused sequence accessions. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
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
print('Removing files with names identical to the output.')
removal = ('rm -r *' + input_ext + 'nice *.faafusednice *' + output_ext + ' ' + sys.argv[4] + ' 2> /dev/null')
os.system(removal)

print('Creating temporary alignment files in single-line FASTA format.')
# modify each multi-line fasta alignment into single-line
for fname2 in filenames:
    #Separate the fasta file stem to use for the temporary single-line alignment files.
    filestem = str(os.path.basename(fname2).split(os.extsep, 1)[0])
    with open(str(datasetsdir + filestem + input_ext), 'r') as multi_input, open(str(filestem + input_ext + 'nice'), 'w') as single_temp:
        block = []
        for line in multi_input:
            if line.startswith('>'):
                if block:
                    single_temp.write(''.join(block) + '\n')
                    block = []
                single_temp.write(line)
            else:
                block.append(line.strip())
        if block:
            single_temp.write(''.join(block) + '\n')

print('Fusing fragmented adjacent sequences and writing accessions to log.')
output_log = open(sys.argv[4], 'w')
for fname3 in filenames:
    filestem = str(os.path.basename(fname3).split(os.extsep, 1)[0])
    output_log.write(filestem + input_ext + '\n')
    outfile = open(str(filestem + '.faafusednice'), 'a')
    # make header-sequence dictionary
    with open(str(filestem + input_ext + 'nice'), 'r') as fasta_file:
        headers = []
        sequences = []
        for line in fasta_file:
            if line.startswith('>'):
                headers.append(line.strip())
            else:
                sequences.append(line.strip())
    headers_sequences = {headers[i]: sequences[i] for i in range(len(headers))}

    # calculate the median length of all sequences
    ungapped_seqs = [s.replace('-','') for s in sequences]
    seqs_lengths = [len(s) for s in ungapped_seqs]
    median_dataset_length = statistics.median(seqs_lengths)

    # make unique contigs list
    with open(str(filestem + input_ext + 'nice'), 'r') as fasta_file:
        unique_contigs = []
        for line in fasta_file:
            match = re.match(r'^>.*? (.*?)_\d+ .*', line)
            if match:
                contig = match.group(1)
                if contig not in unique_contigs:
                    unique_contigs.append(contig)
    # Read the file again to process each unique contig
    with open(str(filestem + input_ext + 'nice'), 'r') as fasta_file:
        for unique_contig in unique_contigs:
            contig_numbers = []
            seq_to_fuse = []
            fasta_file.seek(0)  # Reset file pointer to beginning
            for line in fasta_file:
                match = re.match(r'^>.*? (.*?)_\d+ .*', line)
                if match:
                    if unique_contig == match.group(1):
                        x = line.split(' ')
                        contig_number = int(x[1].split('_')[-1])
                        contig_numbers.append(contig_number)
                        contig_numbers.sort()
            if len(contig_numbers) == 2 and abs(contig_numbers[0] - contig_numbers[1]) == 1:
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[0]))
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[1]))
                leading_gap_count1 = 0
                leading_gap_count2 = 0
                for header, sequence in headers_sequences.items():
                    y = header.split(' ')
                    if y[1] == seq_to_fuse[0]:
                        sequence1 = sequence
                        header1 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count1 += 1
                            else:
                                break
                    elif y[1] == seq_to_fuse[1]:
                        sequence2 = sequence
                        header2 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count2 += 1
                            else:
                                break
                if leading_gap_count1 == leading_gap_count2:
                    outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n')
                else:
                    if leading_gap_count1 > leading_gap_count2:
                        header1_parts = header1.split()
                        header2_parts = header2.split()
                        fused_header = f"{header2_parts[0]} {header2_parts[1]}_{ ' '.join(header1_parts[1:])}"
                        fused_sequence = sequence2 + sequence1
                        fused_sequence = fused_sequence.replace('-', '') # remove all gaps from fused sequence
                    elif leading_gap_count1 < leading_gap_count2:
                        header1_parts = header1.split()
                        header2_parts = header2.split()
                        fused_header = f"{header2_parts[0]} {header1_parts[1]}_{ ' '.join(header2_parts[1:])}"
                        fused_sequence = sequence1 + sequence2
                        fused_sequence = fused_sequence.replace('-', '') # remove all gaps from fused sequence
                    if len(fused_sequence) <= 1.5*median_dataset_length:
                        outfile.write(fused_header + '\n' + fused_sequence + '\n')
                        output_log.write(fused_header + '\n')
                    else:
                        outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n')
            elif len(contig_numbers) == 3 and abs(contig_numbers[0] - contig_numbers[1]) == 1 and abs(contig_numbers[1] - contig_numbers[2]) == 1: # case where all 3 are consecutive
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[0]))
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[1]))
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[2]))
                leading_gap_count1 = 0
                leading_gap_count2 = 0
                leading_gap_count3 = 0
                trailing_gap_count1 = 0
                trailing_gap_count2 = 0
                trailing_gap_count3 = 0
                for header, sequence in headers_sequences.items():
                    y = header.split(' ')
                    if y[1] == seq_to_fuse[0]:
                        sequence1 = sequence
                        header1 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count1 += 1
                            else:
                                break
                        for char in reversed(sequence): # count c-terminal gaps
                            if char == '-':
                                trailing_gap_count1 += 1
                            else:
                                break
                    elif y[1] == seq_to_fuse[1]:
                        sequence2 = sequence
                        header2 = header
                        for char in sequence: # count c-terminal gaps
                            if char == '-':
                                leading_gap_count2 += 1
                            else:
                                break
                        for char in reversed(sequence): # count c-terminal gaps
                            if char == '-':
                                trailing_gap_count2 += 1
                            else:
                                break
                    elif y[1] == seq_to_fuse[2]:
                        sequence3 = sequence
                        header3 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count3 += 1
                            else:
                                break
                        for char in reversed(sequence): # count c-terminal gaps
                            if char == '-':
                                trailing_gap_count3 += 1
                            else:
                                break
                # define first, middle and last sequences for fusion
                if leading_gap_count1 == leading_gap_count2 or leading_gap_count1 == leading_gap_count3 or leading_gap_count2 == leading_gap_count3:
                    outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n' + header3 + '\n' + sequence3 + '\n')
                else:
                    if leading_gap_count1 > leading_gap_count2 and leading_gap_count1 > leading_gap_count3:
                        last_sequence = sequence1
                        last_header = header1
                    elif leading_gap_count2 > leading_gap_count1 and leading_gap_count2 > leading_gap_count3:
                        last_sequence = sequence2
                        last_header = header2
                    elif leading_gap_count3 > leading_gap_count1 and leading_gap_count3 > leading_gap_count2:
                        last_sequence = sequence3
                        last_header = header3
                    if trailing_gap_count1 > trailing_gap_count2 and trailing_gap_count1 > trailing_gap_count3:
                        first_sequence = sequence1
                        first_header = header1
                        if last_sequence == sequence2:
                            middle_sequence = sequence3
                            middle_header = header3
                        else:
                            middle_sequence = sequence2
                            middle_header = header2
                    elif trailing_gap_count2 > trailing_gap_count1 and trailing_gap_count2 > trailing_gap_count3:
                        first_sequence = sequence2
                        first_header = header2
                        if last_sequence == sequence1:
                            middle_sequence = sequence3
                            middle_header = header3
                        else:
                            middle_sequence = sequence1
                            middle_header = header1
                    elif trailing_gap_count3 > trailing_gap_count1 and trailing_gap_count3 > trailing_gap_count2:
                        first_sequence = sequence3
                        first_header = header3
                        if last_sequence == sequence1:
                            middle_sequence = sequence2
                            middle_header = header2
                        else:
                            middle_sequence = sequence1
                            middle_header = header1
                    header1_parts = first_header.split()
                    header2_parts = middle_header.split()
                    header3_parts = last_header.split()
                    fused_header = f"{header2_parts[0]} {header1_parts[1]}_{header2_parts[1]}_{ ' '.join(header3_parts[1:])}"
                    fused_sequence = first_sequence + middle_sequence + last_sequence
                    fused_sequence = fused_sequence.replace('-', '') # remove all gaps from fused sequence
                    if len(fused_sequence) <= 1.5*median_dataset_length:
                        outfile.write(fused_header + '\n' + fused_sequence + '\n')
                        output_log.write(fused_header + '\n')
                    else:
                        outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n' + header3 + '\n' + sequence3 + '\n')
            elif len(contig_numbers) == 3 and abs(contig_numbers[0] - contig_numbers[1]) == 1 and abs(contig_numbers[1] - contig_numbers[2]) != 1: # case where 2 are consecutive (1st-2nd)
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[0]))
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[1]))
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[2])) # 3rd non-consecutive sequence, will not actually get fused
                leading_gap_count1 = 0
                leading_gap_count2 = 0
                for header, sequence in headers_sequences.items():
                    y = header.split(' ')
                    if y[1] == seq_to_fuse[0]:
                        sequence1 = sequence
                        header1 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count1 += 1
                            else:
                                break
                    elif y[1] == seq_to_fuse[1]:
                        sequence2 = sequence
                        header2 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count2 += 1
                            else:
                                break
                    elif y[1] == seq_to_fuse[2]: # write the 3rd non consecutive sequence to the outfile
                        outfile.write(header + '\n' + sequence + '\n')
                if leading_gap_count1 == leading_gap_count2:
                    outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n')
                else:
                    if leading_gap_count1 > leading_gap_count2:
                        header1_parts = header1.split()
                        header2_parts = header2.split()
                        fused_header = f"{header2_parts[0]} {header2_parts[1]}_{ ' '.join(header1_parts[1:])}"
                        fused_sequence = sequence2 + sequence1
                        fused_sequence = fused_sequence.replace('-', '') # remove all gaps from fused sequence
                    elif leading_gap_count1 < leading_gap_count2:
                        header1_parts = header1.split()
                        header2_parts = header2.split()
                        fused_header = f"{header2_parts[0]} {header1_parts[1]}_{ ' '.join(header2_parts[1:])}"
                        fused_sequence = sequence1 + sequence2
                        fused_sequence = fused_sequence.replace('-', '') # remove all gaps from fused sequence
                    if len(fused_sequence) <= 1.5*median_dataset_length:
                        outfile.write(fused_header + '\n' + fused_sequence + '\n')
                        output_log.write(fused_header + '\n')
                    else:
                        outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n')
            elif len(contig_numbers) == 3 and abs(contig_numbers[0] - contig_numbers[1]) != 1 and abs(contig_numbers[1] - contig_numbers[2]) == 1: # case where 2 are consecutive (1st-2nd)
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[1]))
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[2]))
                seq_to_fuse.append(unique_contig + '_' + str(contig_numbers[0])) # 3rd non-consecutive sequence, will not actually get fused
                leading_gap_count1 = 0
                leading_gap_count2 = 0
                for header, sequence in headers_sequences.items():
                    y = header.split(' ')
                    if y[1] == seq_to_fuse[0]:
                        sequence1 = sequence
                        header1 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count1 += 1
                            else:
                                break
                    elif y[1] == seq_to_fuse[1]:
                        sequence2 = sequence
                        header2 = header
                        for char in sequence: # count n-terminal gaps
                            if char == '-':
                                leading_gap_count2 += 1
                            else:
                                break
                    elif y[1] == seq_to_fuse[2]: # write the 3rd non consecutive sequence to the outfile
                        outfile.write(header + '\n' + sequence + '\n')
                if leading_gap_count1 == leading_gap_count2:
                    outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n')
                else:
                    if leading_gap_count1 > leading_gap_count2:
                        header1_parts = header1.split()
                        header2_parts = header2.split()
                        fused_header = f"{header2_parts[0]} {header2_parts[1]}_{ ' '.join(header1_parts[1:])}"
                        fused_sequence = sequence2 + sequence1
                        fused_sequence = fused_sequence.replace('-', '') # remove all gaps from fused sequence
                    elif leading_gap_count1 < leading_gap_count2:
                        header1_parts = header1.split()
                        header2_parts = header2.split()
                        fused_header = f"{header2_parts[0]} {header1_parts[1]}_{ ' '.join(header2_parts[1:])}"
                        fused_sequence = sequence1 + sequence2
                        fused_sequence = fused_sequence.replace('-', '') # remove all gaps from fused sequence
                    if len(fused_sequence) <= 1.5*median_dataset_length:
                        outfile.write(fused_header + '\n' + fused_sequence + '\n')
                        output_log.write(fused_header + '\n')
                        #print(seq_to_fuse)
                    else:
                        outfile.write(header1 + '\n' + sequence1 + '\n' + header2 + '\n' + sequence2 + '\n')
            else:
                for number in contig_numbers: # do a loop for the unlikely case of 4 fragments, better to write all fragments than skip them
                    seq_to_fuse.append(unique_contig + '_' + str(number))
                for item in seq_to_fuse:
                    for header, sequence in headers_sequences.items():
                        y = header.split(' ')
                        if y[1] == item:
                            outfile.write(header + '\n' + sequence + '\n')
    output_log.write('//' + '\n')
    outfile.close()

#Remove gaps from the output FASTA files.
print('Dealigning output FASTA files.')
gap = '-'
for fname4 in filenames:
    filestem = str(os.path.basename(fname4).split(os.extsep, 1)[0])

    with open(str(filestem + output_ext), 'w') as o:
        for record in SeqIO.parse(filestem + '.faafusednice', "fasta"):
            record.seq = record.seq.replace(gap, "")
            SeqIO.write(record, o, "fasta")

print('Removing temporary single-line FASTA files.')
removal2 = ('rm -r *' + input_ext + 'nice *.faafusednice 2> /dev/null')
os.system(removal2)
output_log.close()

print('All done!')
