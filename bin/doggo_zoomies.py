#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This is the WhereDoGGo? wrapper for running phylogenies in IQ-TREE.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)
#2) pandas (https://github.com/pandas-dev/pandas)
#3) NumPy (https://numpy.org/install/)
#4) ETE3 (http://etetoolkit.org/download/ or https://anaconda.org/etetoolkit/ete3)
#5) IQ-TREE 2 (v2.2 or higher) (https://anaconda.org/bioconda/iqtree)

import argparse
import os
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"Bio" : "https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython",
                        "ete3": "https://anaconda.org/etetoolkit/ete3 or https://pypi.org/project/ete3/",
                        "numpy" : "https://numpy.org/install/",
                        "pandas" : "https://github.com/pandas-dev/pandas"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

#Check if required external programs are installed.
import subprocess
externalprograms = {"iqtree2" : "https://anaconda.org/bioconda/iqtree"}
for extprg,link in externalprograms.items():
    try:
        extprg_check = (subprocess.check_output("which " + extprg, shell=True, universal_newlines=True).strip())
    except subprocess.CalledProcessError:
        print('External program ' + extprg + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

#Check if internal scripts are in the PATH.
try:
    fixdualsupports_py = (subprocess.check_output("which fixdualsupports.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script fixdualsupports.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    fixleaves_py = (subprocess.check_output("which fixleaves.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script fixleaves.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    aarecode_py = (subprocess.check_output("which aarecode.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script aarecode.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    desaturate_py = (subprocess.check_output("which desaturate.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script desaturate.py not found in PATH. Exiting.')
    sys.exit(1)
try:
    monophylyspr_py = (subprocess.check_output("which monophylyspr.py", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script monophylyspr.py not found in PATH. Exiting.')
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

parser = argparse.ArgumentParser(description="Henlo, am doggo v20240621. Need halp so I get zoomies?")
parser.add_argument("-i", "--input", required=True, help="INPUT must be the input FASTA file. (required)")
parser.add_argument("-MFP", "--MFP", action='store_true', help="MFP will run IQ-TREE with MODELFINDER to select the model (matrices: LG, WAG, JTT; frequencies: FU,F,FO). (optional)")
parser.add_argument("-C60", "--C60", action='store_true', help="C60 will run IQ-TREE with the C60 mixture model under the matrix picked by MFP (and with the MFP phylogeny as guide tree under the PMSF approximation) or LG (if the MFP log and phylogeny not available) and 10 FreeRate categories i.e., LG (or WAG or JTT )+C60+R10. (optional)")
parser.add_argument("-SR4", "--SR4", action='store_true', help="SR4 will recode data to the 4-state SR alphabet and run IQ-TREE as with the MFP option but with the GTR4 matrix. (optional)")
parser.add_argument("-SR4C60", "--SR4C60", action='store_true', help="SR4C60 will recode data to the 4-state SR alphabet and run IQ-TREE with the SR4C60 model. If available, the SR4 phylogeny will be used as guide tree. (optional)")
parser.add_argument("-D6", "--D6", action='store_true', help="D6 will recode data to the 6-state Dayhoff alphabet and run IQ-TREE as with the MFP option but with the GTR6 matrix. (optional)")
parser.add_argument("-nex", "--nex", required=False, help="NEX must be the input Nexus file for running mixture model analyses on recoded datasets. (optional, unless the SR4C60 option or DESAT option with the SR4C60 argument is selected)")
parser.add_argument("-GHOST", "--GHOST", action='store_true', help="GHOST will run IQ-TREE with the GHOST heterotachy model with the maximum number of mixture categories calculated based on the dataset's number of sequences and positions. (optional)")
parser.add_argument("-desat", "--desat", nargs='+', required=False, help="DESAT will create a set of desaturated datasets at 5%% incements (from 5%% to 90%% of positions retained) for each of the arguments given and run IQ-TREE with the appropriate options. For any argument given, the corresponding option must also be picked or have been run previously to calculate site-specific rates for the desaturation. DESAT arguments can be MFP, C60, SR4, or SR4C60. (optional)")
parser.add_argument("-AU", "--AU", nargs='+', required=False, help="AU will take a user-defined set of leaves and, if they are monophyletic, place them on all possible positions on the tree, optimize branch lengths, and run the Approximately Unbiased test on them versus one or more of the phylogenies run previously. The latter are specified with the same names as their respective options. AU arguments can be MFP, C60, SR4, SR4C60, or D6. (optional)")
parser.add_argument("-AUclade", "--AUclade", required=False, help="AUclade must be a text file with the clade (one leaf name per line) whose position on the phylogeny will be tested. (optional, unless the AU option is provided)")
parser.add_argument("-leaves", "--leaves", required=True, help="LEAVES must be a text file with complete leaf names that will be used to converting them in the phylogenies e.g., an .assembliesnames file. All instances of tab-delimitation will be converted to spaces. (required)")
args=parser.parse_args()

print('Henlo, am doggo v20240621. I get zoomies nao. I speak info messages in hooman lingo.' + '\n')

concat_file_stem = str(os.path.basename(args.input).split(os.extsep, 1)[0])
concat_file_ext = str('.' + str(os.path.basename(args.input).split(os.extsep, 1)[1])) # we have to add the '.' before the extension

# check for input file
check_file = os.path.isfile(args.input)
if check_file == True:
    print('Input file found. Proceeding.')
    #Convert path to an absolute path.
    args.input = os.path.abspath(args.input)
else:
    print('Input file not found. Exiting.')
    sys.exit(1)

# check for leaves file
#TODO: In WhereDoGGo this is expected to be an assembliesnames file. Maybe add a formatting check?
if args.leaves:
    check_file = os.path.isfile(args.leaves)
    if check_file == True:
        print('Input leaf names file found. Proceeding.')
        #Convert path to an absolute path.
        args.leaves = os.path.abspath(args.leaves)
    else:
        print('Input leaf names file not found. Exiting.')
        sys.exit(1)

# check for input nexus file
if args.SR4C60 or (args.desat is not None and "SR4C60" in args.desat):
    check_file = os.path.isfile(args.nex)
    if check_file == True:
        print('Recoded mixture model Nexus file found. Proceeding.')
        #Convert path to an absolute path.
        args.nex = os.path.abspath(args.nex)
    else:
        print('Recoded mixture model Nexus file not found. Exiting.')
        sys.exit(1)

# check for AUclade input file
if args.AU:
    check_file = os.path.isfile(args.AUclade)
    if check_file == True:
        print('Input AUclade file found. Proceeding.')
        #Convert path to an absolute path.
        args.AUclade = os.path.abspath(args.AUclade)
    else:
        print('Input AUclade file not found. Exiting.')
        sys.exit(1)

#Remove any previous output files with the same name.
print ('Removing directories with names identical to the output.')
if args.MFP:
    removal = ('rm -r ' + concat_file_stem + '_MFP/ 2> /dev/null')
    os.system(removal)
if args.C60:
    removal = ('rm -r ' + concat_file_stem + '_C60/ 2> /dev/null')
    os.system(removal)
if args.SR4:
    removal = ('rm -r ' + concat_file_stem + '_SR4/ 2> /dev/null')
    os.system(removal)
if args.D6:
    removal = ('rm -r ' + concat_file_stem + '_D6/ 2> /dev/null')
    os.system(removal)
if args.SR4C60:
    removal = ('rm -r ' + concat_file_stem + '_SR4C60/ 2> /dev/null')
    os.system(removal)
if args.GHOST:
    removal = ('rm -r ' + concat_file_stem + '_GHOST/ 2> /dev/null')
    os.system(removal)
if args.desat is not None:
    for desat_dataset in args.desat:
        removal = ('rm -r ' + concat_file_stem + '_desat_' + desat_dataset + '/ 2> /dev/null')
        os.system(removal)
if args.AU is not None:
    for AU_dataset in args.AU:
        removal = ('rm -r ' + concat_file_stem + '_AU_' + AU_dataset + '/ 2> /dev/null')
        os.system(removal)

#TODO: It might be easier to convert all these commands to for loops.
if args.MFP:
    print('Running phylogeny with the MFP option.')
    run_MFP = str('mkdir ' + concat_file_stem + '_MFP && cp ' + args.input + ' ' + concat_file_stem + '_MFP/ && cd ' + concat_file_stem + '_MFP && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m MFP -mset JTT,WAG,LG -mfreq FU,F,FO -bb 1000 -alrt 1000 -nt AUTO')
    if os.WEXITSTATUS(os.system(run_MFP)) == 1:
        print('Error when running phylogeny with the MFP option. Exiting.')
        sys.exit(1)
    sbs_MFP = str('cd ' + concat_file_stem + '_MFP && python -u ' + fixdualsupports_py + ' ' + concat_file_stem + concat_file_ext + '.treefile ' + concat_file_stem + '_singlesupports.tree >> ' + concat_file_stem + '.fixdualsupportslog')
    if os.WEXITSTATUS(os.system(sbs_MFP)) == 1:
        print('Error when converting branch supports in MFP phylogeny. Exiting.')
        sys.exit(1)
    fl_MFP = str('cd ' + concat_file_stem + '_MFP && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_singlesupports.tree ' + args.leaves + ' ' + concat_file_stem + '_singlesupports_fixedleaves.tree >> ' + concat_file_stem + '.fixleaveslog')
    if os.WEXITSTATUS(os.system(fl_MFP)) == 1:
        print('Error when renaming leaves in MFP phylogeny. Exiting.')
        sys.exit(1)

if args.C60:
    if os.path.isfile(concat_file_stem + '_MFP/' + concat_file_stem + concat_file_ext + '.treefile') and os.path.isfile(concat_file_stem + '_MFP/' + concat_file_stem + concat_file_ext + '.log'):
        print('Running phylogeny with the C60 option with the MFP guide tree.')
        matrix_C60_command = str('grep "Best-fit model" ' + concat_file_stem + '_MFP/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?\\: (.*?)\\+.*/$1/\'')
        matrix_C60 = subprocess.check_output(matrix_C60_command, shell=True, universal_newlines=True).strip()
        optimal_threads_C60_command = str('grep "BEST NUMBER OF THREADS" ' + concat_file_stem + '_MFP/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?: (.*)/$1/\'')
        optimal_threads_C60 = subprocess.check_output(optimal_threads_C60_command, shell=True, universal_newlines=True).strip()
        run_C60 = str('mkdir ' + concat_file_stem + '_C60 && cp ' + args.input + ' ' + concat_file_stem + '_C60/ && cd ' + concat_file_stem + '_C60 && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m ' + matrix_C60 + '+C60+R10 -mwopt -ft ../' + concat_file_stem + '_MFP/' + concat_file_stem + concat_file_ext + '.treefile -bb 1000 -alrt 1000 -nt ' + optimal_threads_C60)
        if os.WEXITSTATUS(os.system(run_C60)) == 1:
            print('Error when running phylogeny with the C60 option. Exiting.')
            sys.exit(1)
        sbs_C60 = str('cd ' + concat_file_stem + '_C60 && python -u ' + fixdualsupports_py + ' ' + concat_file_stem + concat_file_ext + '.treefile ' + concat_file_stem + '_singlesupports.tree >> ' + concat_file_stem + '.fixdualsupportslog')
        if os.WEXITSTATUS(os.system(sbs_C60)) == 1:
            print('Error when converting branch supports in C60 phylogeny. Exiting.')
            sys.exit(1)
        fl_C60 = str('cd ' + concat_file_stem + '_C60 && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_singlesupports.tree ' + args.leaves + ' ' + concat_file_stem + '_singlesupports_fixedleaves.tree >> ' + concat_file_stem + '.fixleaveslog')
        if os.WEXITSTATUS(os.system(fl_C60)) == 1:
            print('Error when renaming leaves in C60 phylogeny. Exiting.')
            sys.exit(1)
    else:
        print('Running phylogeny with the C60 option without a guide tree.')
        run_C60 = str('mkdir ' + concat_file_stem + '_C60 && cp ' + args.input + ' ' + concat_file_stem + '_C60/ && cd ' + concat_file_stem + '_C60 && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m LG+C60+R10 -mwopt -bb 1000 -alrt 1000 -nt AUTO')
        if os.WEXITSTATUS(os.system(run_C60)) == 1:
            print('Error when running phylogeny with the C60 option. Exiting.')
            sys.exit(1)
        sbs_C60 = str('cd ' + concat_file_stem + '_C60 && python -u ' + fixdualsupports_py + ' ' + concat_file_stem + concat_file_ext + '.treefile ' + concat_file_stem + '_singlesupports.tree >> ' + concat_file_stem + '.fixdualsupportslog')
        if os.WEXITSTATUS(os.system(sbs_C60)) == 1:
            print('Error when converting branch supports in C60 phylogeny. Exiting.')
            sys.exit(1)
        fl_C60 = str('cd ' + concat_file_stem + '_C60 && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_singlesupports.tree ' + args.leaves + ' ' + concat_file_stem + '_singlesupports_fixedleaves.tree >> ' + concat_file_stem + '.fixleaveslog')
        if os.WEXITSTATUS(os.system(fl_C60)) == 1:
            print('Error when renaming leaves in C60 phylogeny. Exiting.')
            sys.exit(1)

if args.SR4:
    print('Running phylogeny with the SR4 option.')
    #TODO: Maybe have that file named _SR4.faa but this will cause some necessary rewriting of commands in desat and AU. Ditto for SR4C60 and D6.
    SR4_recode = str('mkdir ' + concat_file_stem + '_SR4 && cd ' + concat_file_stem + '_SR4 && python -u ' + aarecode_py + ' ' + args.input + ' ' + concat_file_stem + concat_file_ext + ' SR4 >> ' + concat_file_stem + '.aarecodelog')
    if os.WEXITSTATUS(os.system(SR4_recode)) == 1:
        print('Error when recoding dataset for SR4 analysis. Exiting.')
        sys.exit(1)
    run_SR4 = str('cd ' + concat_file_stem + '_SR4 && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -mset GTR -mfreq FU,F,FO -bb 1000 -alrt 1000 -nt AUTO')
    if os.WEXITSTATUS(os.system(run_SR4)) == 1:
        print('Error when running phylogeny with the SR4 option. Exiting.')
        sys.exit(1)
    sbs_SR4 = str('cd ' + concat_file_stem + '_SR4 && python -u ' + fixdualsupports_py + ' ' + concat_file_stem + concat_file_ext + '.treefile ' + concat_file_stem + '_singlesupports.tree >> ' + concat_file_stem + '.fixdualsupportslog')
    if os.WEXITSTATUS(os.system(sbs_SR4)) == 1:
        print('Error when converting branch supports in SR4 phylogeny. Exiting.')
        sys.exit(1)
    fl_SR4 = str('cd ' + concat_file_stem + '_SR4 && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_singlesupports.tree ' + args.leaves + ' ' + concat_file_stem + '_singlesupports_fixedleaves.tree >> ' + concat_file_stem + '.fixleaveslog')
    if os.WEXITSTATUS(os.system(fl_SR4)) == 1:
        print('Error when renaming leaves in SR4 phylogeny. Exiting.')
        sys.exit(1)

if args.SR4C60:
    SR4C60_recode = str('mkdir ' + concat_file_stem + '_SR4C60 && cd ' + concat_file_stem + '_SR4C60 && python -u ' + aarecode_py + ' '+ args.input + ' ' + concat_file_stem + concat_file_ext + ' SR4 >> ' + concat_file_stem + '.aarecodelog')
    if os.WEXITSTATUS(os.system(SR4C60_recode)) == 1:
        print('Error when recoding dataset for SR4C60 analysis. Exiting.')
        sys.exit(1)
    if os.path.isfile(concat_file_stem + '_SR4/' + concat_file_stem + concat_file_ext + '.treefile') and os.path.isfile(concat_file_stem + '_SR4/' + concat_file_stem + concat_file_ext + '.log'):
        print('Running phylogeny with the SR4C60 option with the SR4 guide tree.')
        optimal_threads_SR4C60_command = str('grep "BEST NUMBER OF THREADS" ' + concat_file_stem + '_SR4/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?: (.*)/$1/\'')
        optimal_threads_SR4C60 = subprocess.check_output(optimal_threads_SR4C60_command, shell=True, universal_newlines=True).strip()
        run_SR4C60 = str('cd ' + concat_file_stem + '_SR4C60 && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m C60SR4 -mdef ' + args.nex + ' -ft ../'  + concat_file_stem + '_SR4/' + concat_file_stem + concat_file_ext + '.treefile -bb 1000 -alrt 1000 -nt ' + optimal_threads_SR4C60)
        if os.WEXITSTATUS(os.system(run_SR4C60)) == 1:
            print('Error when running phylogeny with the SR4C60 option. Exiting.')
            sys.exit(1)
        sbs_SR4C60 = str('cd ' + concat_file_stem + '_SR4C60 && python -u ' + fixdualsupports_py + ' ' + concat_file_stem + concat_file_ext + '.treefile ' + concat_file_stem + '_singlesupports.tree >> ' + concat_file_stem + '.fixdualsupportslog')
        if os.WEXITSTATUS(os.system(sbs_SR4C60)) == 1:
            print('Error when converting branch supports in SR4C60 phylogeny. Exiting.')
            sys.exit(1)
        fl_SR4C60 = str('cd ' + concat_file_stem + '_SR4C60 && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_singlesupports.tree ' + args.leaves + ' ' + concat_file_stem + '_singlesupports_fixedleaves.tree >> ' + concat_file_stem + '.fixleaveslog')
        if os.WEXITSTATUS(os.system(fl_SR4C60)) == 1:
            print('Error when renaming leaves in SR4C60 phylogeny. Exiting.')
            sys.exit(1)
    else:
        print('Running phylogeny with the SR4C60 option without a guide tree.')
        run_SR4C60 = str('cd ' + concat_file_stem + '_SR4C60 && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m C60SR4 -mdef ' + args.nex + ' -bb 1000 -alrt 1000 -nt AUTO')
        if os.WEXITSTATUS(os.system(run_SR4C60)) == 1:
            print('Error when running phylogeny with the SR4C60 option. Exiting.')
            sys.exit(1)
        sbs_SR4C60 = str('cd ' + concat_file_stem + '_SR4C60 && python -u ' + fixdualsupports_py + ' ' + concat_file_stem + concat_file_ext + '.treefile ' + concat_file_stem + '_singlesupports.tree >> ' + concat_file_stem + '.fixdualsupportslog')
        if os.WEXITSTATUS(os.system(sbs_SR4C60)) == 1:
            print('Error when converting branch supports in SR4C60 phylogeny. Exiting.')
            sys.exit(1)
        fl_SR4C60 = str('cd ' + concat_file_stem + '_SR4C60 && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_singlesupports.tree ' + args.leaves + ' ' + concat_file_stem + '_singlesupports_fixedleaves.tree >> ' + concat_file_stem + '.fixleaveslog')
        if os.WEXITSTATUS(os.system(fl_SR4C60)) == 1:
            print('Error when renaming leaves in SR4C60 phylogeny. Exiting.')
            sys.exit(1)

if args.D6:
    print('Running phylogeny with the D6 option.')
    D6_recode = str('mkdir ' + concat_file_stem + '_D6 && cd ' + concat_file_stem + '_D6 && python -u ' + aarecode_py + ' ' + args.input + ' ' + concat_file_stem + concat_file_ext + ' D6 >> ' + concat_file_stem + '.aarecodelog')
    if os.WEXITSTATUS(os.system(D6_recode)) == 1:
        print('Error when recoding dataset for D6 analysis. Exiting.')
        sys.exit(1)
    run_D6 = str('cd ' + concat_file_stem + '_D6 && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -mset GTR -mfreq FU,F,FO -bb 1000 -alrt 1000 -nt AUTO')
    if os.WEXITSTATUS(os.system(run_D6)) == 1:
        print('Error when running phylogeny with the D6 option. Exiting.')
        sys.exit(1)
    sbs_D6 = str('cd ' + concat_file_stem + '_D6 && python -u ' + fixdualsupports_py + ' ' + concat_file_stem + concat_file_ext + '.treefile ' + concat_file_stem + '_singlesupports.tree >> ' + concat_file_stem + '.fixdualsupportslog')
    if os.WEXITSTATUS(os.system(sbs_D6)) == 1:
        print('Error when converting branch supports in D6 phylogeny. Exiting.')
        sys.exit(1)
    fl_D6 = str('cd ' + concat_file_stem + '_D6 && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_singlesupports.tree ' + args.leaves + ' ' + concat_file_stem + '_singlesupports_fixedleaves.tree >> ' + concat_file_stem + '.fixleaveslog')
    if os.WEXITSTATUS(os.system(fl_D6)) == 1:
        print('Error when renaming leaves in D6 phylogeny. Exiting.')
        sys.exit(1)

if args.GHOST:
    print('Running phylogeny with the GHOST option. WARNING: In this option no singlesupports and fixedleaves trees will be produced.')
    # This is for cmax calculation
    A = 0
    n = 0
    #Get number of positions
    with open(args.input, "r") as alignment:
        for line in alignment:
            if line.startswith(">"):  # Skip header lines if present
                continue
            A = len(line.strip())  # Count characters in the sequence
            n += 1
    # Calculate the maximum value of k (cmax)
    k = int(((A/10)+1) // ((2 * n) + 17)) ## the // returns the min integer
    cmax = k
    cmax = str(cmax)
    print(cmax)
    run_GHOST = str('mkdir ' + concat_file_stem + '_GHOST && cp ' + args.input + ' ' + concat_file_stem + '_GHOST/ && cd ' + concat_file_stem + '_GHOST && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -mset JTT,WAG,LG -mrate H,*H -cmax ' + cmax + ' -mfreq F,FU,FO -bb 1000 -alrt 1000 -nt AUTO')
    if os.WEXITSTATUS(os.system(run_GHOST)) == 1:
        print('Error when running phylogeny with the GHOST option. Exiting.')
        sys.exit(1)

if args.desat is not None:
    desat_options = ['MFP', 'C60', 'SR4', 'SR4C60']
    for desat_dataset in args.desat:
        if desat_dataset in desat_options:
            pass
        else:
            print('Desaturation dataset ' + desat_dataset + ' not valid. Exiting.')
            sys.exit(1)
    print('All desaturation datasets are valid. Proceeding.')
    for desat_dataset in args.desat:
        if os.path.isfile(concat_file_stem + '_' + desat_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile') and os.path.isfile(concat_file_stem + '_' + desat_dataset + '/' + concat_file_stem + concat_file_ext + '.log'):
            pass
        else:
            print('Corresponding phylogeny or log file for ' + desat_dataset + ' not found. Exiting.')
            sys.exit(1)
    print('All corresponding phylogenies and log files for desaturation datasets found. Proceeding.')
    for desat_dataset in args.desat:
        print('Running desaturation analysis for the ' + desat_dataset + ' option.')
        setup_directories = str('mkdir ' + concat_file_stem + '_desat_' + desat_dataset + ' ' + concat_file_stem + '_desat_' + desat_dataset + '/rates_Poisson_G16 && cp ' + args.input + ' ' + concat_file_stem + '_desat_' + desat_dataset + '/rates_Poisson_G16/')
        if os.WEXITSTATUS(os.system(setup_directories)) == 1:
            print('Error when setting up directories for ' + desat_dataset + ' desaturation analysis. Exiting.')
            sys.exit(1)
        get_rates = str('cd ' + concat_file_stem + '_desat_'  + desat_dataset + '/rates_Poisson_G16 && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m Poisson+G16 -n 0 -t ../../' +  concat_file_stem + '_' + desat_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile -blfix --rate -nt AUTO')
        if os.WEXITSTATUS(os.system(get_rates)) == 1:
            print('Error when calculating site-specific rates for ' + desat_dataset + ' desaturation analysis. Exiting.')
            sys.exit(1)
        create_datasets = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && python -u ' + desaturate_py + ' rates_Poisson_G16/' + concat_file_stem + concat_file_ext + '.rate 20 rates_Poisson_G16/' + concat_file_stem + concat_file_ext + ' >> ' + concat_file_stem + '.desaturatelog')
        if os.WEXITSTATUS(os.system(create_datasets)) == 1:
            print('Error when creating desaturated datasets for ' + desat_dataset + ' analysis. Exiting.')
            sys.exit(1)
        if desat_dataset == "SR4" or desat_dataset == "SR4C60":
            recode_datasets = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && mkdir original_datasets && mv *.faa original_datasets/ && for i in original_datasets/*.faa ; do python -u ' + aarecode_py + ' $i "$(basename $i .faa)".faa SR4 >> "$(basename $i .faa)".aarecodelog ; done')
            if os.WEXITSTATUS(os.system(recode_datasets)) == 1:
                print('Error when recoding desaturated datasets for ' + desat_dataset + ' analysis. Exiting.')
                sys.exit(1)
        if desat_dataset== "MFP":
            run_desaturation = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && for i in *.faa ; do iqtree2 -quiet -s $i -m MFP -mset JTT,WAG,LG -mfreq FU,F,FO -bb 1000 -alrt 1000 -nt AUTO ; done')
        if desat_dataset== "C60":
            matrix_C60_command = str('grep "Command\\:" ' + concat_file_stem + '_C60/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*? -m (.*?)\\+.*/$1/\'')
            matrix_C60 = subprocess.check_output(matrix_C60_command, shell=True, universal_newlines=True).strip()
            check_threads_C60_command = str('grep "Command\\:" ' + concat_file_stem + '_C60/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*? -nt (.*)/$1/\'')
            check_threads_C60 = subprocess.check_output(check_threads_C60_command, shell=True, universal_newlines=True).strip()
            if check_threads_C60 == "AUTO":
                optimal_threads_C60_command = str('grep "BEST NUMBER OF THREADS" ' + concat_file_stem + '_C60/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?: (.*)/$1/\'')
                optimal_threads_C60 = subprocess.check_output(optimal_threads_C60_command, shell=True, universal_newlines=True).strip()
            else:
                optimal_threads_C60 = check_threads_C60
            run_desaturation = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && for i in *.faa ; do iqtree2 -quiet -s $i -m ' + matrix_C60 + '+C60+R10 -mwopt -ft ../' + concat_file_stem + '_C60/' + concat_file_stem + concat_file_ext + '.treefile -bb 1000 -alrt 1000 -nt ' + optimal_threads_C60 + ' ; done')
        if desat_dataset== "SR4":
            run_desaturation = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && for i in *.faa ; do iqtree2 -quiet -s $i -mset GTR -mfreq FU,F,FO -bb 1000 -alrt 1000 -nt AUTO ; done')
        if desat_dataset== "SR4C60":
            check_threads_SR4C60_command = str('grep "Command\\:" ' + concat_file_stem + '_SR4C60/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*? -nt (.*)/$1/\'')
            check_threads_SR4C60 = subprocess.check_output(check_threads_SR4C60_command, shell=True, universal_newlines=True).strip()
            if check_threads_SR4C60 == "AUTO":
                optimal_threads_SR4C60_command = str('grep "BEST NUMBER OF THREADS" ' + concat_file_stem + '_SR4C60/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?: (.*)/$1/\'')
                optimal_threads_SR4C60 = subprocess.check_output(optimal_threads_SR4C60_command, shell=True, universal_newlines=True).strip()
            else:
                optimal_threads_SR4C60 = check_threads_SR4C60
            run_desaturation = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && for i in *.faa ; do iqtree2 -quiet -s $i -m C60SR4 -mdef ' + args.nex + ' -ft ../' + concat_file_stem + '_SR4C60/' + concat_file_stem + concat_file_ext + '.treefile -bb 1000 -alrt 1000 -nt ' + optimal_threads_SR4C60 + ' ; done')
        if os.WEXITSTATUS(os.system(run_desaturation)) == 1:
            print('Error when running phylogenies with the ' + desat_dataset + ' option on desaturated datasets. Exiting.')
            sys.exit(1)
        sbs_desaturation = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && for i in *.faa.treefile ; do python -u ' + fixdualsupports_py + ' $i "$(basename $i .faa.treefile)"_singlesupports.tree >> "$(basename $i .faa.treefile)".fixdualsupportslog ; done')
        if os.WEXITSTATUS(os.system(sbs_desaturation)) == 1:
            print('Error when converting branch supports in desaturated dataset ' + desat_dataset + ' phylogenies. Exiting.')
            sys.exit(1)
        fl_desaturation = str('cd ' + concat_file_stem + '_desat_' + desat_dataset + ' && for i in *_singlesupports.tree ; do python -u ' + fixleaves_py + ' $i ' + args.leaves + ' "$(basename $i .tree)"_fixedleaves.tree >> "$(basename $i .tree | perl -p -e \'s/^(.*?)_singlesupports/$1/\')".fixleaveslog ; done')
        if os.WEXITSTATUS(os.system(fl_desaturation)) == 1:
            print('Error when renaming leaves in desaturated dataset ' + desat_dataset + ' phylogenies. Exiting.')
            sys.exit(1)

if args.AU is not None:
    AU_options = ['MFP', 'C60', 'SR4', 'SR4C60', 'D6']
    for AU_dataset in args.AU:
        if AU_dataset in AU_options:
            pass
        else:
            print('AU dataset ' + AU_dataset + ' not valid. Exiting.')
            sys.exit(1)
    print('All AU datasets are valid. Proceeding.')
    for AU_dataset in args.AU:
        if os.path.isfile(concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile') and os.path.isfile(concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + '.log'):
            pass
        else:
            print('Corresponding phylogeny or log file for ' + AU_dataset + ' not found. Exiting.')
            sys.exit(1)
    print('All corresponding phylogenies and log files for AU datasets found. Proceeding.')
    for AU_dataset in args.AU:
        setup_directories = str('mkdir ' + concat_file_stem + '_AU_' + AU_dataset + ' && cp ' + concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + ' ' + concat_file_stem + '_AU_' + AU_dataset + '/')
        if os.WEXITSTATUS(os.system(setup_directories)) == 1:
            print('Error when setting up directories for ' + AU_dataset + ' AU test. Exiting.')
            sys.exit(1)
        print('Running AU topology test for the ' + AU_dataset + ' option.')
        run_monophylyspr = str('cd ' + concat_file_stem + '_AU_' + AU_dataset + ' && python -u ' + monophylyspr_py + ' ../'  + concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + '_singlesupports.tree ' + args.AUclade + ' >> ' + concat_file_stem + '.monophylysprlog')
        if os.WEXITSTATUS(os.system(run_monophylyspr)) == 1:
            print('Error when running monophylyspr.py for ' + AU_dataset + ' AU test. Exiting.')
            sys.exit(1)
        if AU_dataset == "MFP":
            model_MFP_command = str('grep "Best-fit model" ' + concat_file_stem + '_MFP/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?\\: (.*?) .*/$1/\'')
            model_MFP = subprocess.check_output(model_MFP_command, shell=True, universal_newlines=True).strip()
            run_AUtest = str('cd ' + concat_file_stem + '_AU_' + AU_dataset + ' && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m ' + model_MFP + ' -zb 10000 -z SPR_trees.treels -te ../' + concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile -au -nt AUTO')
        if AU_dataset == "C60":
            matrix_C60_command = str('grep "Command\\:" ' + concat_file_stem + '_C60/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*? -m (.*?)\\+.*/$1/\'')
            matrix_C60 = subprocess.check_output(matrix_C60_command, shell=True, universal_newlines=True).strip()
            run_AUtest = str('cd ' + concat_file_stem + '_AU_' + AU_dataset + ' && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m ' + matrix_C60 + '+C60+R10 -mwopt -zb 10000 -z SPR_trees.treels -te ../' + concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile -au -nt AUTO')
        if AU_dataset == "SR4":
            model_SR4_command = str('grep "Best-fit model" ' + concat_file_stem + '_SR4/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?\\: (.*?) .*/$1/\'')
            model_SR4 = subprocess.check_output(model_SR4_command, shell=True, universal_newlines=True).strip()
            run_AUtest = str('cd ' + concat_file_stem + '_AU_' + AU_dataset + ' && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m ' + model_SR4 + ' -zb 10000 -z SPR_trees.treels -te ../' + concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile -au -nt AUTO')
        if AU_dataset == "SR4C60":
            run_AUtest = str('cd ' + concat_file_stem + '_AU_' + AU_dataset + ' && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m C60SR4 -mdef ' + args.nex + ' -zb 10000 -z SPR_trees.treels -te ../' + concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile -au -nt AUTO')
        if AU_dataset == "D6":
            model_D6_command = str('grep "Best-fit model" ' + concat_file_stem + '_D6/' + concat_file_stem + concat_file_ext + '.log | perl -p -e \'s/^.*?\\: (.*?) .*/$1/\'')
            model_D6 = subprocess.check_output(model_D6_command, shell=True, universal_newlines=True).strip()
            run_AUtest = str('cd ' + concat_file_stem + '_AU_' + AU_dataset + ' && iqtree2 -quiet -s ' + concat_file_stem + concat_file_ext + ' -m ' + model_D6 + ' -zb 10000 -z SPR_trees.treels -te ../' + concat_file_stem + '_' + AU_dataset + '/' + concat_file_stem + concat_file_ext + '.treefile -au -nt AUTO')
        if os.WEXITSTATUS(os.system(run_AUtest)) == 1:
            print('Error when running AU test for the ' + AU_dataset + ' dataset. Exiting.')
            sys.exit(1)
        fl_AUtest = str('cd ' + concat_file_stem + '_AU_' + AU_dataset + ' && perl -p -e \'s/^\\[.*?\\](.*)/$1/g\' ' + concat_file_stem + concat_file_ext + '.trees >> ' + concat_file_stem + '_nolikelihoods.trees && python -u ' + fixleaves_py + ' ' + concat_file_stem + '_nolikelihoods.trees ' + args.leaves + ' ' + concat_file_stem + '_nolikelihoods_fixedleaves_nosupports.trees >> ' + concat_file_stem + '.fixleaveslog && perl -p -i -e \'s/\\)\\d+\\:/\\)\\:/g\' ' + concat_file_stem + '_nolikelihoods_fixedleaves_nosupports.trees')
        if os.WEXITSTATUS(os.system(fl_AUtest)) == 1:
            print('Error when renaming leaves in AU dataset ' + AU_dataset + ' phylogenies. Exiting.')
            sys.exit(1)

print('Bork bork! I finish. Boop the snoot pls?')
