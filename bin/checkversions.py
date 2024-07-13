#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This checks the version of all WhereDoGGo? dependencies.

#Dependencies
#1) Biopython (https://biopython.org/wiki/Download or https://anaconda.org/conda-forge/biopython)
#2) pandas (https://github.com/pandas-dev/pandas)
#3) NumPy (https://numpy.org/install/)
#4) ETE3 (http://etetoolkit.org/download/ or https://anaconda.org/etetoolkit/ete3)
#5) ncbi-datasets-cli (https://github.com/ncbi/datasets)
#6) Pyrodigal (https://github.com/althonos/pyrodigal)
#7) HMMER (https://anaconda.org/bioconda/hmmer)
#8) Pullseq (https://anaconda.org/bioconda/pullseq)
#9) MAFFT (https://anaconda.org/bioconda/mafft)
#10) BMGE (https://anaconda.org/bioconda/bmge)
#11) IQ-TREE 2 (v2.2 or higher) (https://anaconda.org/bioconda/iqtree)

#NOTE 1: All code was written and tested on Intel macOS with some testing on Linux. Please report any issues.

import importlib
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
externalprograms = {"datasets":"https://github.com/ncbi/datasets",
                    "pyrodigal": "https://github.com/althonos/pyrodigal",
                    "hmmsearch": "https://anaconda.org/bioconda/hmmer",
                    "pullseq": "https://anaconda.org/bioconda/pullseq",
                    "einsi": "https://anaconda.org/bioconda/mafft",
                    "bmge": "https://anaconda.org/bioconda/bmge",
                    "iqtree2" : "https://anaconda.org/bioconda/iqtree"}
for extprg,link in externalprograms.items():
    try:
        extprg_check = (subprocess.check_output("which " + extprg, shell=True, universal_newlines=True).strip())
    except subprocess.CalledProcessError:
        print('External program ' + extprg + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

print('#Script: checkversions.py')
print('#Version: v20240713')
print('#Usage: python checkversions.py')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 1:
    print ('No arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

print(f"Python {sys.version.split()[0]}")

nonstandardlibraries_list = ["biopython", "ete3", "numpy", "pandas"] #We have to take this indirect route, because pip lists the library as biopython, not as Bio
for nstl_v in nonstandardlibraries_list:
    nonstandardlibraries_version = (subprocess.check_output("pip list --format=freeze | grep " + nstl_v + " | perl -p -e \'s/^.*?==(.*?)/$1/\'", shell=True, universal_newlines=True).strip())
    print(nstl_v + ' ' + nonstandardlibraries_version)

datasets_check = (subprocess.check_output("datasets --version | perl -p -e \'s/^.*?\\: (.*?)/$1/\'", shell=True, universal_newlines=True).strip())
print('datasets ' + datasets_check)

pyrodigal_check = (subprocess.check_output("pyrodigal -V", shell=True, universal_newlines=True).strip())
print(pyrodigal_check)

hmmsearch_check = (subprocess.check_output("hmmsearch -h | grep \"# HMMER\" | perl -p -e \'s/^.*? HMMER (.*?) .*/$1/\'", shell=True, universal_newlines=True).strip())
print('hmmsearch ' + hmmsearch_check)

pullseq_check = (subprocess.check_output("pullseq -i goat --version 2> /dev/null | perl -p -e \'s/^Version is (.*)/$1/\'", shell=True, universal_newlines=True).strip())
print('pullseq ' + pullseq_check)

einsi_check = (subprocess.check_output("einsi --version 2>&1 >/dev/null | perl -p -e \'s/^v(.*?) .*/$1/\'", shell=True, universal_newlines=True).strip())
print('einsi ' + einsi_check)

bmge_check = (subprocess.check_output("bmge -? | grep \"BMGE\" | grep \"version\" | perl -p -e \'s/BMGE \\(version (.*?)\\).*/$1/\'", shell=True, universal_newlines=True).strip())
print('bmge ' + bmge_check)

iqtree2_check = (subprocess.check_output("iqtree2 --version | grep \"IQ-TREE \" | perl -p -e \'s/^IQ-TREE .*? version (.*?) .*/$1/\'", shell=True, universal_newlines=True).strip())
print('iqtree2 ' + iqtree2_check)

print ('WhereDoGGo? modules each use the following dependencies: fetch (datasets, pyrodigal), herd (pyrodigal), sniff (biopython, hmmsearch, pullseq, einsi, bmge), zoomies (biopython, ete3, numpy, pandas, iqtree2).')

print('All done!')
