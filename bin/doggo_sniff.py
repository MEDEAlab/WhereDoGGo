#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This is the WhereDoGGo? wrapper for creating concatenated marker datasets from local genome databases (homology searches, parsing, removing multiplicates, aligning, trimming, concatenation).

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)
#2) HMMER (https://anaconda.org/bioconda/hmmer)
#3) seqtk (https://anaconda.org/bioconda/seqtk)
#4) MAFFT (https://anaconda.org/bioconda/mafft)
#5) BMGE (https://anaconda.org/bioconda/bmge)

import argparse
import os
import random
import re
import string
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

#Check if required external programs are installed.
import subprocess
externalprograms = {"hmmsearch": "https://anaconda.org/bioconda/hmmer",
                    "seqtk": "https://anaconda.org/bioconda/seqtk",
                    "einsi": "https://anaconda.org/bioconda/mafft",
                    "bmge": "https://anaconda.org/bioconda/bmge"}
for extprg,link in externalprograms.items():
    try:
        extprg_check = (subprocess.check_output("which " + extprg, shell=True, universal_newlines=True).strip())
    except subprocess.CalledProcessError:
        print('External program ' + extprg + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

#Check if internal scripts are in the PATH.
try:
    hmmsearchout2accessions_sh = (subprocess.check_output("which hmmsearchout2accessions.sh", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script hmmsearchout2accessions.sh not found in PATH. Exiting.')
    sys.exit(1)
try:
    fasta2distribution_py = (subprocess.check_output("which fasta2distribution.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script fasta2distribution.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    fuseadjacent_py = (subprocess.check_output("which fuseadjacent.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script fuseadjacent.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    removemultiples_py = (subprocess.check_output("which removemultiples.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script removemultiples.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    preconcatenation_sh = (subprocess.check_output("which preconcatenation.sh", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script preconcatenation.sh not found in PATH. Exiting.')
    sys.exit(1)
try:
    concatenation_py = (subprocess.check_output("which concatenation.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script concatenation.py not found in PATH. Exiting.')
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

parser = argparse.ArgumentParser(description="Henlo, am doggo v20241212. Need halp for sniff markers?")
parser.add_argument("-db", "--databases", nargs='+', required=True, help="DATABASES must be one or more multi-FASTA files with genome amino acid sequences (e.g., from the output of doggo_fetch or doggo_herd), against which the HMM searches will be run. (required)")
parser.add_argument("-hmm", "--hmm", required=True, help="HMM must be a directory containing the .hmm files (HMM profiles) for the markers. (trailing slash optional) (required)")
parser.add_argument("-con", "--concatenation", required=False, help="CONCATENATION must be the filename stem of the concatenation output containing alphanumeric characters and/or underscores only. If not provided, it will default to five random alphanumeric characters. (optional)")
parser.add_argument("-cut", "--cutoffs", required=False, help="CUTOFFS must be a tab-delimited file with two columns, the marker HMM profile filename and its domain bitscore cutoff as it would be input in HMMER. The HMM profile filename must only contain alphanumeric characters and/or underscores and use the .hmm extension. The domain bitscore must be a number, with or without decimals. For any files not included or if this argument is not provided, default domain bitscore cutoff is 30. (optional)")
parser.add_argument("-f", "--fuse", action='store_true', help="FUSE will run an additional step to fuse fragmented adjacent sequences (up to three fragments, four or more will be ignored), based on their accessions. WARNING: This option is still experimental, use it with caution and manually check your final alignments and concatenation. (optional)")
args=parser.parse_args()
#TODO: Add possibility for the user to define an output directory.

print('Henlo, am doggo v20241212. I sniff markers nao. I speak info messages in hooman lingo.' + '\n')

## checkpoint for input db file(s)
for dbmember in args.databases:
    #dbmember = os.path.abspath(dbmember)
    if os.path.isfile(dbmember) == True:
        pass
    else:
        print('Database file ' + dbmember + ' not found. Exiting.')
        sys.exit(1)
print('All database files found. Proceeding.')

## checkpoint for hmm directory existence and trailing slash
if os.path.exists(args.hmm) == True:
    print ('Directory with HMM profiles found. Proceeding.')
    args.hmm = os.path.abspath(args.hmm)
    args.hmm = os.path.join(args.hmm, '')
else:
    print ('Directory with HMM profiles not found. Exiting.')
    sys.exit(1)

## checkpoint if hmm directory contains .hmm files
#TODO: Check if those are real hmm profiles by looking at their formatting?
for fname in os.listdir(args.hmm):
    if fname.endswith('.hmm'):
        print ('File(s) with .hmm extension found in the HMM profiles directory. Proceeding.')
        break
else:
    print('No files with the .hmm extension found in the HMM profiles directory. Exiting.')
    sys.exit(1)

## checkpoint if CUTOFFS file exists
if args.cutoffs is not None:
    if os.path.isfile(args.cutoffs) == True:
        print('Cutoffs file found. Proceeding.')
        args.cutoffs = os.path.abspath(args.cutoffs)
    else:
        print('Cutoffs file not found. Exiting.')
        sys.exit(1)

## checkpoint if CUTOFFS file has correct formatting
    with open(args.cutoffs) as f:
        for i, line in enumerate(f):
            lineContains = line.split('\t')
            lineLength = len(lineContains)
            if lineLength != 2 or re.match(r"^[A-Za-z0-9_]+[.]hmm$", lineContains[0]) == False or re.match(r"^\d+[.]?\d*$", lineContains[1]) == False:
                print('Wrongly formatted line ' + line + ' found in cutoffs file. Exiting.')
                sys.exit(1)
    print('Cutoffs file correctly formatted. Proceeding.')
else:
    print('No cutoffs file specified. All HMM searches will use the default domain bitscore cutoff 30. Proceeding.')

#Checkpoint for concatenation name.
if args.concatenation is None:
    print('No concatenation name provided, so one will be randomly generated. Proceeding.')
    args.concatenation = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(5))
elif re.match(r'^[A-Za-z0-9_]+$', args.concatenation):
    print ('Concatenation name is valid. Proceeding.')
else:
    print ('Concatenation name is invalid. Exiting.')
    sys.exit(1)

#Remove any previous output files with the same name.
print ('Removing files and directories with names identical to the output.')
removal = str('rm -r ' + args.concatenation + '_hmmsearch/ ' + args.concatenation + '_hmmsearchout2accessions/ ' + args.concatenation + '_seqtk/ ' + args.concatenation + '_faafixedheaders/ ' + args.concatenation + '_einsiprefuse/ ' + args.concatenation + '_fuseadjacent/ ' + args.concatenation + '_removemultiples/ ' + args.concatenation + '_einsi/ ' + args.concatenation + '_bmge30/ ' + args.concatenation + '_preconcatenation/ ' + args.concatenation + '_otherlogs/ ' + args.concatenation + '_hmmsearch.tar.gz ' + args.concatenation + '_hmmsearchout2accessions.tar.gz ' + args.concatenation + '_seqtk.tar.gz ' + args.concatenation + '_faafixedheaders.tar.gz ' + args.concatenation + '_einsiprefuse.tar.gz ' + args.concatenation + '_fuseadjacent.tar.gz ' + args.concatenation + '_removemultiples.tar.gz ' + args.concatenation + '_einsi.tar.gz ' + args.concatenation + '_bmge30.tar.gz ' + args.concatenation + '_preconcatenation.tar.gz ' + args.concatenation + '_otherlogs.tar.gz ' + args.concatenation + '.database ' + args.concatenation + '.assembliesnames *.distro ' + args.concatenation + '.distribution ' + args.concatenation + '.fasta2distributionlog ' + args.concatenation + '.concatenation ' + args.concatenation + '.concatenationlog ' + args.concatenation + '_sniff/  2> /dev/null')
os.system(removal)

##Create a combined database in the current directory to avoid spreading datasets over multiple directories.
#TODO: Have the hmmsearch run separately for each db and output in a new directory. Then the pipeline runs separately for the output of each db and the fasta files are combined at some point, possibly at demultiplied, before aligning.
print ('Combining individual databases and assembliesnames files.')
for dbmember in args.databases:
    cat_databases = str('cat ' + dbmember + ' >> ' + args.concatenation + '.database')
    if os.WEXITSTATUS(os.system(cat_databases)) == 1:
        print('Error while combining databases. Exiting.')
        sys.exit(1)
    #TODO: Parameter expansion might not be the safest (most portable) way to remove the extensions. Maybe look into dirname combined and basename.
    #cat_assembliesnames = str('cat "${' + dbmember + '%.*}".assembliesnames >> ' + args.concatenation + '.assembliesnames')
    cat_assembliesnames = str('cat "$(dirname ' + dbmember + ')"/"$(basename ' + dbmember + ' | perl -p -e \'s/^(.*?)\\..*/$1/g\')".assembliesnames >> ' + args.concatenation + '.assembliesnames')
    if os.WEXITSTATUS(os.system(cat_assembliesnames)) == 1:
        print('WARNING: No corresponding assembliesnames file found in the same directory and with the same filename stem as ' + dbmember + '.')

#TODO: Maybe all these steps should be a Shell script. It looks weird to have perl one-liners in a Python wrapper.
print ('Running HMM searches.')
set_cutoff = '30'
dirhmmsearch = str('mkdir ' + args.concatenation + '_hmmsearch')
if os.WEXITSTATUS(os.system(dirhmmsearch)) == 1:
    print('Error when making hmmsearch directory.')
    sys.exit(1)
if args.cutoffs is not None:
    correspond_cutoff = {}
    with open(args.cutoffs, 'r') as cutoffs:
        for line in cutoffs:
            x = line.split('\t')
            correspond_cutoff.update({ x[0] : x[1].strip() })
    for fname in os.listdir(args.hmm):
        if fname.endswith('.hmm') and fname in correspond_cutoff.keys():
            set_cutoff = correspond_cutoff[fname]
            hmmsearch=str('cd ' + args.concatenation + '_hmmsearch && hmmsearch --domT ' + set_cutoff + ' --domtblout "$(basename ' + fname + ' .hmm)".hmmsearchout ' + args.hmm + fname + ' ../' + args.concatenation + '.database >> ' + args.concatenation + '.hmmsearchlog 2>&1')
            if os.WEXITSTATUS(os.system(hmmsearch)) == 1:
                print('Error during hmmsearch for ' + args.hmm + fname + '. Exiting.')
                sys.exit(1)
        elif fname.endswith('.hmm') and fname not in correspond_cutoff.keys():
            set_cutoff = '30'
            hmmsearch=str('cd ' + args.concatenation + '_hmmsearch && hmmsearch --domT ' + set_cutoff + ' --domtblout "$(basename ' + fname + ' .hmm)".hmmsearchout ' + args.hmm + fname + ' ../' + args.concatenation + '.database >> ' + args.concatenation + '.hmmsearchlog 2>&1')
            if os.WEXITSTATUS(os.system(hmmsearch)) == 1:
                print('Error during hmmsearch for ' + args.hmm + fname + '. Exiting.')
                sys.exit(1)
else:
    for fname in os.listdir(args.hmm):
        if fname.endswith('.hmm'):
            hmmsearch=str('cd ' + args.concatenation + '_hmmsearch && hmmsearch --domT ' + set_cutoff + ' --domtblout "$(basename ' + fname + ' .hmm)".hmmsearchout ' + args.hmm + fname + ' ../' + args.concatenation + '.database >> ' + args.concatenation + '.hmmsearchlog 2>&1')
            if os.WEXITSTATUS(os.system(hmmsearch)) == 1:
                print('Error during hmmsearch for ' + args.hmm + fname + '. Exiting.')
                sys.exit(1)

## this is for extracting the accessions from the .hmmsearchout (hmm search output)
print ('Extracting marker accessions from HMM search output.')
exacc = str('mkdir ' + args.concatenation + '_hmmsearchout2accessions && cd ' + args.concatenation + '_hmmsearchout2accessions && bash ' + hmmsearchout2accessions_sh + ' ../' + args.concatenation + '_hmmsearch/ .hmmsearchout >> ' + args.concatenation + '.hmmsearchout2accessionslog')
if os.WEXITSTATUS(os.system(exacc)) == 1:
    print('Error during hmmsearchout2accessions.sh script. Exiting.')
    sys.exit(1)

## this is now for pulling the actual sequences from the local database(s)
print ('Pulling sequences from database (seqtk).')
#TODO: seqtk doesn't seem to output anything to stdout or stderror, so there's nothing we can redirect as log.
seqtk = str('mkdir ' + args.concatenation + '_seqtk && cd ' + args.concatenation + '_seqtk && for i in ../' + args.concatenation + '_hmmsearchout2accessions/*.accessions ; do seqtk subseq ../' + args.concatenation + '.database "$i" >> "$(basename $i .accessions)".faaoriginal ; done')
if os.WEXITSTATUS(os.system(seqtk)) == 1:
    print('Error when pulling sequences from the local database. Exiting.')
    sys.exit(1)

## this is for creating a file with the taxonomic distribution of each marker
print ('Creating taxonomic distribution file.')
dirotherlogs = str('mkdir ' + args.concatenation + '_otherlogs')
if os.WEXITSTATUS(os.system(dirotherlogs)) == 1:
    print('Error when making the otherlogs directory. Exiting.')
    sys.exit(1)
taxdistro = str ('python -u ' + fasta2distribution_py + ' ./' + args.concatenation + '_seqtk/ .faaoriginal ' + args.concatenation + '.distribution ' + args.concatenation + '.assembliesnames >> ' + args.concatenation + '.fasta2distributionlog && mv ' + args.concatenation + '.fasta2distributionlog ' + args.concatenation + '_otherlogs/')
if os.WEXITSTATUS(os.system(taxdistro)) == 1:
    print('Error during fasta2distribution.py script. Exiting.')
    sys.exit(1)

## this is for removing weird symbols from the .faaoriginal fasta files that often cause alignment or tree-building programs to crash and reversing the assembly and sequence accessions.
#The reason is to have the same assembly accession as FASTA header accession for all sequences from the same genome across markers, in order to be able to concatenate them.
#There's nothing to redirect as a log here.
print ('Fixing FASTA headers (removing problematic characters, reversing assembly and sequence accessions).')
#Originally was using cp with the new extensions and in-place editing, but after testing, runtime is 1/3 with this command, plus slightly lower peak memory consumption.
headerfix = str('mkdir ' + args.concatenation + '_faafixedheaders && cd ' + args.concatenation + '_faafixedheaders && for i in ../' + args.concatenation + '_seqtk/*.faaoriginal ; do perl -p -e \'s/\\(/ /g\' "$i" | perl -p -e \'s/\\)/ /g\' | perl -p -e \'s/\\;/ /g\' | perl -p -e \'s/\\:/ /g\' | perl -p -e \'s/\\,/ /g\' | perl -p -e \'s/>(.*?) (.*?) (.*)/>$2 $1 $3/g\' >> "$(basename $i .faaoriginal)".faafixedheaders ; done')
if os.WEXITSTATUS(os.system(headerfix)) == 1:
    print('Error when fixing FASTA headers. Exiting.')
    sys.exit(1)

if args.fuse:
    ## this is for fusing adjacent fragmented sequences
    print ('Aligning with MAFFT E-INS-i and fusing adjacent fragmented sequences.')
    #The MAFFT screen output is actually stderror.
    einsiunfused = str('mkdir ' + args.concatenation + '_einsiprefuse && cd ' + args.concatenation + '_einsiprefuse && for i in ../' + args.concatenation + '_faafixedheaders/*.faafixedheaders ; do echo "$(basename $i .faafixedheaders)".faafixedheaders >> ' +  args.concatenation + '.einsiunfusedlog ; einsi --thread -1 --reorder $i > "$(basename $i .faafixedheaders)".einsiunfused 2>> ' +  args.concatenation + '.einsiunfusedlog ; echo "//" >> '  +  args.concatenation + '.einsiunfusedlog ; done')
    if os.WEXITSTATUS(os.system(einsiunfused)) == 1:
        print('Error when aligning datasets with MAFFT E-INS-i (before fusing adjacent fragmented sequences). Exiting.')
        sys.exit(1)
    fuseadjacent = str('mkdir ' + args.concatenation + '_fuseadjacent && cd ' + args.concatenation + '_fuseadjacent && python -u ' + fuseadjacent_py + ' ../' + args.concatenation + '_einsiprefuse/ .einsiunfused .faafused ' + args.concatenation + '.fusedlog >> '  + args.concatenation + '.fuseadjacentlog')
    if os.WEXITSTATUS(os.system(fuseadjacent)) == 1:
        print('Error when fusing adjacent fragmented sequences. Exiting.')
        sys.exit(1)
    ## this is for removing taxa with multiple sequences from datasets
    print ('Removing from each marker any taxa with multiple sequences.')
    removemultiples = str('mkdir ' + args.concatenation + '_removemultiples && cd ' + args.concatenation + '_removemultiples && python -u ' + removemultiples_py + ' ../' + args.concatenation + '_fuseadjacent/ .faafused .faademultiplied ' + args.concatenation + '.demultipliedlog >> '  + args.concatenation + '.removemultipleslog')
else:
    print ('Removing from each marker any taxa with multiple sequences.')
    removemultiples = str('mkdir ' + args.concatenation + '_removemultiples && cd ' + args.concatenation + '_removemultiples && python -u ' + removemultiples_py + ' ../' + args.concatenation + '_faafixedheaders/ .faafixedheaders .faademultiplied ' + args.concatenation + '.demultipliedlog >> '  + args.concatenation + '.removemultipleslog')
if os.WEXITSTATUS(os.system(removemultiples)) == 1:
    print('Error during removemultiples.py script. Exiting.')
    sys.exit(1)

## this is for aligning the .demultiplied fasta files with mafft
print ('Aligning with MAFFT E-INS-i.')
#The MAFFT screen output is actually stderror.
einsi = str('mkdir ' + args.concatenation + '_einsi && cd ' + args.concatenation + '_einsi && for i in ../' + args.concatenation + '_removemultiples/*.faademultiplied ; do echo "$(basename $i .faademultiplied)".faademultiplied >> ' +  args.concatenation + '.einsilog ; einsi --thread -1 --reorder $i > "$(basename $i .faademultiplied)".einsi 2>> ' +  args.concatenation + '.einsilog ; echo "//" >> '  +  args.concatenation + '.einsilog ; done')
if os.WEXITSTATUS(os.system(einsi)) == 1:
    print('Error when aligning datasets with MAFFT E-INS-i. Exiting.')
    sys.exit(1)

## this is for trimming with BMGE
print ('Trimming with BMGE (BLOSUM30).')
bmge30 = ('mkdir ' + args.concatenation + '_bmge30 && cd ' + args.concatenation + '_bmge30 && for i in ../' + args.concatenation + '_einsi/*.einsi ; do bmge -i $i -t AA -m BLOSUM30 -of "$(basename $i .einsi)".bmge30 | perl -p -e \'s/\\r//g\' | perl -p -e \'s/^.*?problem/problem/g\' | perl -p -e \'s/^.*?Amino/Amino/g\' | perl -p -e \'s/^.*?before/before/g\' | perl -p -e \'s/^.*?after/after/g\' | perl -p -e \'s/\\s+after.*//g\' >> ' + args.concatenation + '.bmge30log ; done')
if os.WEXITSTATUS(os.system(bmge30)) == 1:
    print('Error when trimming alignments with BMGE. Exiting.')
    sys.exit(1)

## this is for running the preconcatenation script (responsible for generating a file with the dataset names to be concatenated in the next step)
print ('Running preconcatenation script.')
preconcatenation = str('mkdir ' + args.concatenation + '_preconcatenation && cd ' + args.concatenation + '_preconcatenation && bash ' + preconcatenation_sh + ' ../' + args.concatenation + '_bmge30/ .bmge30 >> ' + args.concatenation + '.preconcatenationlog')
if os.WEXITSTATUS(os.system(preconcatenation)) == 1:
    print('Error during preconcatenation.sh script. Exiting.')
    sys.exit(1)

## this is for running the concatenation script and creating the final .concatenation fasta file with concatenated marker genes
print ('Concatenating markers.')
concatenation = str('python -u ' + concatenation_py + ' ./' + args.concatenation + '_bmge30/ ./' + args.concatenation + '_preconcatenation/dataset.pass ' + args.concatenation + '.concatenation >> ' + args.concatenation + '.concatenationlog && mv ' + args.concatenation + '.concatenationlog ' + args.concatenation + '_otherlogs/')
if os.WEXITSTATUS(os.system(concatenation)) == 1:
    print('Error during concatenation.py script. Exiting.')
    sys.exit(1)

#Back up files in a dedicated directory. Remove the combined database to avoid redundancy and save disk space.
print ('Creating run directory and removing combined database.')
if args.fuse:
    backup = str('mkdir ' + args.concatenation + '_sniff && tar -czf ' + args.concatenation + '_hmmsearch.tar.gz ' + args.concatenation + '_hmmsearch/ && tar -czf ' + args.concatenation +  '_hmmsearchout2accessions.tar.gz ' + args.concatenation + '_hmmsearchout2accessions/ && tar -czf ' + args.concatenation + '_seqtk.tar.gz ' + args.concatenation + '_seqtk/ && tar -czf ' + args.concatenation + '_faafixedheaders.tar.gz ' + args.concatenation + '_faafixedheaders/ && tar -czf ' + args.concatenation + '_einsiprefuse.tar.gz ' + args.concatenation + '_einsiprefuse/ && tar -czf ' + args.concatenation + '_fuseadjacent.tar.gz ' + args.concatenation + '_fuseadjacent/ && tar -czf ' + args.concatenation + '_removemultiples.tar.gz ' + args.concatenation + '_removemultiples/ && tar -czf ' + args.concatenation + '_einsi.tar.gz ' + args.concatenation + '_einsi/ && tar -czf ' + args.concatenation + '_bmge30.tar.gz ' + args.concatenation + '_bmge30/ && tar -czf ' + args.concatenation + '_preconcatenation.tar.gz ' + args.concatenation + '_preconcatenation/ && tar -czf ' + args.concatenation + '_otherlogs.tar.gz ' + args.concatenation + '_otherlogs/ && mv -i ' + args.concatenation + '_hmmsearch.tar.gz ' + args.concatenation + '_hmmsearchout2accessions.tar.gz ' + args.concatenation + '_seqtk.tar.gz ' + args.concatenation + '_faafixedheaders.tar.gz ' + args.concatenation + '_einsiprefuse.tar.gz ' + args.concatenation + '_fuseadjacent.tar.gz ' + args.concatenation + '_removemultiples.tar.gz ' + args.concatenation + '_einsi.tar.gz ' + args.concatenation + '_bmge30.tar.gz ' + args.concatenation + '_preconcatenation.tar.gz ' + args.concatenation + '_otherlogs.tar.gz ' + args.concatenation + '.assembliesnames ' + args.concatenation + '.distribution ' + args.concatenation + '.concatenation ' + args.concatenation + '_sniff/ && rm -r ' + args.concatenation + '_hmmsearch/ ' + args.concatenation + '_hmmsearchout2accessions/ ' + args.concatenation + '_seqtk/ ' + args.concatenation + '_faafixedheaders/ ' + args.concatenation + '_einsiprefuse/ ' + args.concatenation + '_fuseadjacent/ ' + args.concatenation + '_removemultiples/ ' + args.concatenation + '_einsi/ ' + args.concatenation + '_bmge30/ ' + args.concatenation + '_preconcatenation/ ' + args.concatenation + '_otherlogs/ ' + args.concatenation + '.database 2> /dev/null')
else:
    backup = str('mkdir ' + args.concatenation + '_sniff && tar -czf ' + args.concatenation + '_hmmsearch.tar.gz ' + args.concatenation + '_hmmsearch/ && tar -czf ' + args.concatenation + '_hmmsearchout2accessions.tar.gz ' + args.concatenation + '_hmmsearchout2accessions/ && tar -czf ' + args.concatenation + '_seqtk.tar.gz ' + args.concatenation + '_seqtk/ && tar -czf ' + args.concatenation + '_faafixedheaders.tar.gz ' + args.concatenation + '_faafixedheaders/ && tar -czf ' + args.concatenation + '_removemultiples.tar.gz ' + args.concatenation + '_removemultiples/ && tar -czf ' + args.concatenation + '_einsi.tar.gz ' + args.concatenation + '_einsi/ && tar -czf ' + args.concatenation + '_bmge30.tar.gz ' + args.concatenation + '_bmge30/ && tar -czf ' + args.concatenation + '_preconcatenation.tar.gz ' + args.concatenation + '_preconcatenation/ && tar -czf ' + args.concatenation + '_otherlogs.tar.gz ' + args.concatenation + '_otherlogs/ && mv -i ' + args.concatenation + '_hmmsearch.tar.gz ' + args.concatenation + '_hmmsearchout2accessions.tar.gz ' + args.concatenation + '_seqtk.tar.gz ' + args.concatenation + '_faafixedheaders.tar.gz ' + args.concatenation + '_removemultiples.tar.gz ' + args.concatenation + '_einsi.tar.gz ' + args.concatenation + '_bmge30.tar.gz ' + args.concatenation + '_preconcatenation.tar.gz ' + args.concatenation + '_otherlogs.tar.gz ' + args.concatenation + '.assembliesnames ' + args.concatenation + '.distribution ' + args.concatenation + '.concatenation ' + args.concatenation + '_sniff/ && rm -r ' + args.concatenation + '_hmmsearch/ ' + args.concatenation + '_hmmsearchout2accessions/ ' + args.concatenation + '_seqtk/ ' + args.concatenation + '_faafixedheaders/ ' + args.concatenation + '_removemultiples/ ' + args.concatenation + '_einsi/ ' + args.concatenation + '_bmge30/ ' + args.concatenation + '_preconcatenation/ ' + args.concatenation + '_otherlogs/ ' + args.concatenation + '.database 2> /dev/null')
if os.WEXITSTATUS(os.system(backup)) == 1:
    print('Error when creating run directory and removing combined database. Exiting.')
    sys.exit(1)

print('Bork bork! I finish. Gib treato pls?')
