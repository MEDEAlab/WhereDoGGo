#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This is the WhereDoGGo? wrapper for picking a set of genomes, downloading them from Genbank as contigs, and creating a local genomic database of protein sequences.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#1) ncbi-datasets-cli (https://github.com/ncbi/datasets)
#2) Pyrodigal (https://github.com/althonos/pyrodigal)

import argparse
import math
import os
import sys

#Check if required external programs are installed.
import subprocess
externalprograms = {"datasets":"https://github.com/ncbi/datasets",
                    "pyrodigal": "https://github.com/althonos/pyrodigal"}
for extprg,link in externalprograms.items():
    try:
        extprg_check = (subprocess.check_output("which " + extprg, shell=True, universal_newlines=True).strip())
    except subprocess.CalledProcessError:
        print('External program ' + extprg + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

#Check if internal scripts are in the PATH.
try:
    pickgenomes_py = (subprocess.check_output("which pickgenomes.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script pickgenomes.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    downloadcontigs_sh = (subprocess.check_output("which downloadcontigs.sh", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script downloadcontigs.sh not found in PATH. Exiting.')
    sys.exit(1)
try:
    contigs2orfs_sh = (subprocess.check_output("which contigs2orfs.sh", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script contigs2orfs.sh not found in PATH. Exiting.')
    sys.exit(1)
try:
    createdb_sh = (subprocess.check_output("which createdb.sh", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script createdb.sh not found in PATH. Exiting.')
    sys.exit(1)

print(r"""
                           ▄              ▄
                          ▌▒█           ▄▀▒▌
                          ▌▒▒█        ▄▀▒▒▒▐
                         ▐▄▀▒▒▀▀▀▀▄▄▄▀▒▒▒▒▒▐
                       ▄▄▀▒░▒▒▒▒▒▒▒▒▒█▒▒▄█▒▐
                     ▄▀▒▒▒░░░▒▒▒░░░▒▒▒▀██▀▒▌
                    ▐▒▒▒▄▄▒▒▒▒░░░▒▒▒▒▒▒▒▀▄▒▒▌
                    ▌░░▌█▀▒▒▒▒▒▄▀█▄▒▒▒▒▒▒▒█▒▐
                   ▐░░░▒▒▒▒▒▒▒▒▌██▀▒▒░░░▒▒▒▀▄▌
                   ▌░▒▄██▄▒▒▒▒▒▒▒▒▒░░░░░░▒▒▒▒▌
                  ▌▒▀▐▄█▄█▌▄░▀▒▒░░░░░░░░░░▒▒▒▐
                  ▐▒▒▐▀▐▀▒░▄▄▒▄▒▒▒▒▒▒░▒░▒░▒▒▒▒▌
                  ▐▒▒▒▀▀▄▄▒▒▒▄▒▒▒▒▒▒▒▒░▒░▒░▒▒▐
                   ▌▒▒▒▒▒▒▀▀▀▒▒▒▒▒▒░▒░▒░▒░▒▒▒▌
                   ▐▒▒▒▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▒▄▒▒▐
                    ▀▄▒▒▒▒▒▒▒▒▒▒▒░▒░▒░▒▄▒▒▒▒▌
                      ▀▄▒▒▒▒▒▒▒▒▒▒▄▄▄▀▒▒▒▒▄▀
                        ▀▄▄▄▄▄▄▀▀▀▒▒▒▒▒▄▄▀
                           ▒▒▒▒▒▒▒▒▒▒▀▀
__        ___                   ____         ____  ____      ___
\ \      / / |__   ___ _ __ ___|  _ \  ___  / ___|/ ___| ___|__ \
 \ \ /\ / /| '_ \ / _ \ '__/ _ \ | | |/ _ \| |  _| |  _ / _ \ / /
  \ V  V / | | | |  __/ | |  __/ |_| | (_) | |_| | |_| | (_) |_|
   \_/\_/  |_| |_|\___|_|  \___|____/ \___/ \____|\____|\___/(_)
                    """)

parser = argparse.ArgumentParser(description="Henlo, am doggo v20240713. Need halp for fetch genomes?")
parser.add_argument("-i", "--input", required=True, help="INPUT must be a tab-delimited file as per GTDB\'s metadata files for Bacteria or Archaea. (required)")
parser.add_argument("-lvl", "--level", required=True, help="LEVEL must be the highest taxonomic level for genomes to be selected as presented in GTDB taxonomy strings e.g., p__Asgardarchaeota (i.e., get genomes from within LEVEL). (required)")
parser.add_argument("-res", "--resolution", required=True, help="RESOLUTION must be a taxonomic level lower than LEVEL. For each RESOLUTION in LEVEL, NUMBER genomes are picked (or as many as available). Must be p (phylum), c (class), o (order), f (family), g (genus), or all. (required)")
parser.add_argument("-n", "--number", required=True, help="NUMBER must be the number of genomes to be picked per RESOLUTION. Must be a positive integer or all. If RESOLUTION is all, NUMBER must also be all. (required)")
parser.add_argument("-ig", "--ignore", required=False, help="IGNORE must be a text file containing genome assembly accessions (one per line, versionless e.g., GCA_011362025) that will not bepicked. (optional)")
args=parser.parse_args()
file_stem = (str(args.level)+'_'+str(args.resolution)+'_'+str(args.number))

print('Henlo, am doggo v20240713. I fetch genomes nao. I speak info messages in hooman lingo.' + '\n')

# check for input file
#TODO: Include formatting check for INPUT and IGNORE.
if os.path.isfile(args.input) == True:
    print('Input file found. Proceeding.')
else:
    print('Input file not found. Exiting.')
    sys.exit(1)

#Check if input file is parsed GTDB metadata.
with open(args.input, 'r') as parsedcheck:
    for line in parsedcheck:
        x = line.split('\t')
        if x[18]=='t':
            pass
        else:
            print('Input file does not contain parsed GTDB metadata. Exiting.')
            sys.exit(1)
print('Input file contains parsed GTDB metadata. Proceeding.')

# check for taxonomic level
with open(args.input, 'r') as file:
    content = file.read()
    if str(args.level+';') in content:
        print('Taxonomic level is valid. Proceeding.')
    else:
        print('Taxonomic level is invalid. Exiting.')
        sys.exit(1)

# checkpoint for resolution and number of genomes both being "all"
if args.resolution == 'all' and args.number != 'all':
    print('When taxonomic resolution is "all", number must also be "all". Exiting.')
    sys.exit(1)
elif args.resolution != 'all' and args.number == 'all':
    print('When number is "all", taxonomic resolution must also be "all". Exiting.')
    sys.exit(1)

# check for taxonomic resolution
tax_res = ['p', 'c', 'o', 'f', 'g', 'all']
if args.resolution in tax_res:
    print('Taxonomic resolution is valid. Proceeding.')
else:
    print('Taxonomic resolution is invalid. Exiting.')
    sys.exit(1)

# check for number of genomes to be picked
try:
    args.number=int(args.number)
except:
    pass
if (isinstance(args.number, int) and int(args.number)>0) or (args.number == 'all'):
    print('Number given is valid. Proceeding.')
else:
    print('Number given is invalid. Exiting.')
    sys.exit(1)
args.number=str(args.number)

# check for ignore_list argument
if args.ignore is not None:
    #check_ignore = os.path.isfile(path+'/'+args.ignore)
    if os.path.isfile(args.ignore) == True:
        print('Ignore list file found. Proceeding.')
    else:
        print('Ignore list file not found. Exiting.')
        sys.exit(1)
else:
        print('No ignore list specified. Proceeding.')

#Remove any previous output files with the same name.
print ('Removing files and directories with names identical to the output.')
removal = str('rm -r ' + file_stem + '.assemblies ' + file_stem + '.assembliesnames ' + file_stem + '.pickgenomeslog ' + file_stem + '.downloadcontigslog ' + file_stem + '.failed ' + file_stem + '.contigs2orfslog ' + file_stem + '.pyrodigallog ' + file_stem + '.createdblog ' + file_stem + '.database ' + file_stem + '_contigs.tar.gz ' + file_stem + '_orfs.tar.gz ' + file_stem + '_contigs/ ' + file_stem + '_orfs/ ' + file_stem + '_fetch/ 2> /dev/null')
os.system(removal)

# this is for defining the domain of the taxonomic level chosen by the user
with open(args.input, 'r') as f:
    for line in f:
        x = line.split('\t')
        z = x[19].split(';')
        if args.level in x[19]:
            if z[0] == 'd__Archaea':
                domain = 'Archaea'
                break
            elif z[0] == 'd__Bacteria':
                domain = 'Bacteria'
                break

print ('Picking genomes.')
if len(sys.argv) == 11:
    pickgenomes = str('python -u ' + pickgenomes_py + ' ' + args.input + ' ' + args.level + ' ' + args.resolution + ' ' + args.number + ' ' + args.ignore + ' >> ' + file_stem + '.pickgenomeslog')
    #os.system(pickgenomes)
    if os.WEXITSTATUS(os.system(pickgenomes)) == 1:
        print('Error during pickgenomes.py script. Exiting.')
        sys.exit(1)
else:
    pickgenomes = str('python -u ' + pickgenomes_py + ' ' + args.input + ' ' + args.level + ' ' + args.resolution + ' ' + args.number+ ' >> ' + file_stem + '.pickgenomeslog')
    #os.system(pickgenomes)
    if os.WEXITSTATUS(os.system(pickgenomes)) == 1:
        print('Error during pickgenomes.py script. Exiting.')
        sys.exit(1)

print ('Downloading contigs. Additional download rounds are retries for corrupted files.')
download = str('bash ' + downloadcontigs_sh + ' ' + file_stem + '.assemblies >> ' + file_stem + '.downloadcontigslog')
#os.system(download)
if os.WEXITSTATUS(os.system(download)) in (1, 5):
    print('Error or warning during downloadcontigs.sh script. Exiting.')
    sys.exit(1)

print ('Predicting ORFs from contigs.')
contigs2orfs = str('bash ' + contigs2orfs_sh + ' ' + file_stem + '.assemblies ' + file_stem + '_contigs/ ' + '.fna assemblies >> ' + file_stem + '.contigs2orfslog')
#os.system(contigs2orfs)
if os.WEXITSTATUS(os.system(contigs2orfs)) == 1:
    print('Error during contigs2orfs.sh script. Exiting.')
    sys.exit(1)

print ('Creating local database.')
createdb = str('bash ' + createdb_sh + ' ' + file_stem + '.assemblies >> ' + file_stem + '.createdblog')
#os.system(createdb)
if os.WEXITSTATUS(os.system(createdb)) == 1:
    print('Error during createdb.sh script. Exiting.')
    sys.exit(1)

print ('Compressing contigs and ORFs directories.')
compressdirs = str('tar -czf ' + file_stem + '_contigs.tar.gz ' + file_stem + '_contigs/ && rm -r ' + file_stem + '_contigs/ && tar -czf ' + file_stem + '_orfs.tar.gz ' + file_stem + '_orfs/ && rm -r ' + file_stem + '_orfs/ 2> /dev/null')
if os.WEXITSTATUS(os.system(compressdirs)) == 1:
    print('Error when compressing contigs and ORFs directories. Exiting.')
    sys.exit(1)

print ('Creating run directory.')
backup = str('mkdir ' + file_stem + '_fetch && mv -i ' + file_stem + '.assemblies ' + file_stem + '.assembliesnames ' + file_stem + '.pickgenomeslog ' + file_stem + '.downloadcontigslog ' + file_stem + '.failed ' + file_stem + '.contigs2orfslog ' + file_stem + '.pyrodigallog ' + file_stem + '.createdblog ' + file_stem + '.database ' + file_stem + '_contigs.tar.gz ' + file_stem + '_orfs.tar.gz ' + file_stem + '_fetch/' )
if os.WEXITSTATUS(os.system(backup)) == 1:
    print('Error when creating run directory. Exiting.')
    sys.exit(1)

print('Bork bork! I finish. Gib chimken pls?')
