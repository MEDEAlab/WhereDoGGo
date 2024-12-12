#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This is the WhereDoGGo? wrapper for creating a local genomic protein sequence database from a set of genome files (contigs).

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.

#Dependencies
#1) Pyrodigal (https://github.com/althonos/pyrodigal)

import argparse
import os
import random
import re
import string
import sys

#Check if required external programs are installed.
import subprocess
externalprograms = {"pyrodigal":"https://github.com/althonos/pyrodigal"}
for extprg,link in externalprograms.items():
    try:
        extprg_check = (subprocess.check_output("which " + extprg, shell=True, universal_newlines=True).strip())
    except subprocess.CalledProcessError:
        print('External program ' + extprg + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

try:
    createrecords_sh = (subprocess.check_output("which createrecords.sh", shell=True, universal_newlines=True).strip())
except subprocess.CalledProcessError:
    print('Script createrecords.sh not found in PATH. Exiting.')
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

parser = argparse.ArgumentParser(description="Henlo, am doggo v20241212. Need halp for local genomes?")
parser.add_argument("-loc", "--location", required=True, help='LOCATION must be a directory containing all genome files. The trailing slash will be added, if not included. (required)')
parser.add_argument("-ext", "--extension", required=True, help='EXTENSION must be the file extension of the genome files. The leading dot will be added, if not included. (required)')
parser.add_argument("-prj", "--project", required=False, help='PROJECT must be a project code consisting of alphanumeric characters and/or underscores only. This code will be used as the filename stem of the generated files. If not provided, it will default to five random alphanumeric characters. (optional)')
parser.add_argument("-code25", "--code25", action='store_true', help="CODE25 must be provided if the local genomes are SR1 or Gracilibacteria (c__JAEDAM01). All ORF predictions will use Genetic Code 25. (optional)")
args=parser.parse_args()

print('Henlo, am doggo v20241212. I prepare local genomes nao. I speak info messages in hooman lingo.' + '\n')

# checkpoint for path existence and trailing slash
if os.path.exists(args.location) == True:
    print ('Genome directory found. Proceeding.')
    args.location = os.path.join(args.location, '')
else:
    print ('Genome directory not found. Exiting.')
    sys.exit(1)

# checkpoint for file extension argument leading dot
if args.extension.startswith('.') == False:
    args.extension = str('.' + args.extension)

# checkpoint if files with given extension exist
for fname in os.listdir(args.location):
    if fname.endswith(args.extension):
        print ('File(s) with the given extension found in the genome directory. Proceeding.')
        break
else:
    print('No files with the given extension found in the genome directory. Exiting.')
    sys.exit(1)

# checkpoint project name
if args.project is None:
    print('No project code given, so one will be randomly generated. Proceeding.')
    args.project = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(5))
elif re.match(r'^[A-Za-z0-9_]+$', args.project):
    print('Project code is valid. Proceeding.')
else:
    print('Project code is invalid. Exiting.')
    sys.exit(1)

#Remove any previous output files with the same name.
print ('Removing files and directories with names identical to the output.')
removal = ('rm -r ' + args.project + '.assemblies ' + args.project + '.assembliesnames ' + args.project + '.createrecordslog ' + args.project + '.contigs2orfslog ' + args.project + '.pyrodigallog ' + args.project + '.createdblog ' + args.project + '.database ' + args.project + '_orfs.tar.gz' + args.project + '_orfs/ ' + args.project + '_herd/ 2> /dev/null')
os.system(removal)

print ('Creating records for local genomes.')
createrecords = str('bash ' + createrecords_sh + ' ' + args.location + ' ' + args.extension + ' ' + args.project + ' >> '  + args.project + '.createrecordslog')
#os.system(createrecords)
if os.WEXITSTATUS(os.system(createrecords)) == 1:
    print('Error during createrecords.sh script. Exiting.')
    sys.exit(1)

print ('Predicting ORFs from contigs.')
if args.code25:
    contigs2orfs = str('bash ' + contigs2orfs_sh + ' ' + args.project + '.assemblies ' + args.location + ' ' + args.extension + ' names code25 >> '  + args.project + '.contigs2orfslog')
else:
    contigs2orfs = str('bash ' + contigs2orfs_sh + ' ' + args.project + '.assemblies ' + args.location + ' ' + args.extension + ' names >> '  + args.project + '.contigs2orfslog')
#os.system(contigs2orfs)
if os.WEXITSTATUS(os.system(contigs2orfs)) == 1:
    print('Error during contigs2orfs.sh script. Exiting.')
    sys.exit(1)

print ('Creating local database.')
createdb = str('bash ' + createdb_sh + ' ' + args.project + '.assemblies '  + args.project + '_orfs/ ' + '.faa >> ' + args.project + '.createdblog')
#os.system(createdb)
if os.WEXITSTATUS(os.system(createdb)) == 1:
    print('Error during createdb.sh script. Exiting.')
    sys.exit(1)

print ('Compressing ORFs directory.')
compressdirs = str('tar -czf ' + args.project + '_orfs.tar.gz ' + args.project + '_orfs/ && rm -r ' + args.project + '_orfs/ 2> /dev/null')
if os.WEXITSTATUS(os.system(compressdirs)) == 1:
    print('Error when compressing contigs and ORFs directories. Exiting.')
    sys.exit(1)

print ('Creating run directory.')
backup = str('mkdir ' + args.project + '_herd && mv -i ' + args.project + '.assemblies ' + args.project + '.assembliesnames ' + args.project + '.createrecordslog ' + args.project + '.contigs2orfslog ' + args.project + '.pyrodigallog ' + args.project + '.createdblog ' + args.project + '.database ' + args.project + '_orfs.tar.gz ' + args.project + '_herd/')
if os.WEXITSTATUS(os.system(backup)) == 1:
    print('Error when creating run directory. Exiting.')
    sys.exit(1)

print('Bork bork! I finish. Gib chimken pls?')
