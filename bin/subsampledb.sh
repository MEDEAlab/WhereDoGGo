#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script subsamples from a local database a set of genomes (e.g., picked with pickgenomes or from a custom list). Both assemblies found and not found are logged. It's meant mainly as a utility script for subsampling, if one keeps around complete domain-level databases, to save time from running doggo_fetch.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.
#NOTE 2: This is the WhereDoGGo? version of the script that assumes input files are correctly formatted, as per the output of doggo_fetch.

#Dependencies
#1) seqtk (https://anaconda.org/bioconda/seqtk)

assemblies="$1"
inputdb="$2"

#Check if all dependencies are installed.
if ! command -v seqtk /dev/null
then
    echo "Program seqtk not installed. Download it from https://anaconda.org/bioconda/seqtk. Exiting."
    exit 1
fi

cat << EndOfMessage
#Script: subsampledb.sh
#Version: v20241212
#Usage: subsampledb.sh <assemblies> <inputdb>
#<assemblies> must be a text file of versionless assemblies that will be subsampled from inputdb (1/line). (required)
#<inputdb> must be a local database in FASTA format. (required)
#For more information refer to the comments in the script and/or the Github page.
EndOfMessage

#Check if the number of arguments is correct, otherwise exit.
if [ "$#" -eq 2 ]
then
	echo "Two arguments found. Proceeding."
else
	echo "Wrong number of arguments given. Exiting."
	exit 1
fi

#Check if the assemblies file exists, otherwise exit.
if [ -f $assemblies ]
then
	echo "Assemblies file found. Proceeding."
else
	echo "Assemblies file not found. Exiting."
	exit 1
fi

#Check if the inputdb file exists, otherwise exit.
if [ -f $inputdb ]
then
	echo "Local database file found. Proceeding."
else
	echo "Local database file not found. Exiting."
	exit 1
fi

#Isolate the stem of the assemblies file name as a separate variable to avoid any issues with extensions.
baseassemblies="$(basename "$assemblies" | perl -p -e 's/^(.*?)\..*/$1/g')"

#For any given assemblies file, remove all the files a previous run would have produced to avoid clashes.
echo "Removing files and directories with names identical to the output."
rm -r "$baseassemblies".database "$baseassemblies".subsampled "$baseassemblies".notsubsampled 2> /dev/null

#Assuming the FASTA headers are in doggo format (>accession assembly [species]), reverse the assembly and accession, so that each assembly can be used to pull all sequences from the same genome.
echo "Reversing accessions and assemblies in the local database FASTA headers. WARNING: If the script crashes now, the database will be unusable. Refer to the manual on how to fix this issue."
perl -p -i -e 's/>(.*?) (.*?) (.*)/>$2 $1 $3/g' "$inputdb"

#Pull the sequences from inputdb.
echo "Pulling sequences from the local database."
seqtk subseq "$inputdb" "$assemblies" >> "$baseassemblies".database

#Check which assemblies were subsampled.
echo "Checking which assemblies were subsampled."
perl -n -e 'print if m/^>/g' "$baseassemblies".database | perl -p -e 's/>(.*?) .*/$1/g' | uniq >> "$baseassemblies".subsampled

#Reverse the assembly and accession in inputdb and the subsampled db to the original order.
echo "Reversing accessions and assemblies in the local and subsampled database FASTA headers. WARNING: If the script crashes now, the databases will be unusable. Refer to the manual on how to fix this issue."
perl -p -i -e 's/>(.*?) (.*?) (.*)/>$2 $1 $3/g' "$inputdb"
perl -p -i -e 's/>(.*?) (.*?) (.*)/>$2 $1 $3/g' "$baseassemblies".database

#Check which of the assemblies were not in the inputdb (not subsampled).
echo "Checking which assemblies were not subsampled."
while IFS= read -r line ; do if LC_ALL=C fgrep -q -- "$line" "$baseassemblies".subsampled ; then : ; else echo "$line" >> "$baseassemblies".notsubsampled ; fi ; done < "$assemblies"

#Congrats, you're done!
echo "All done!"
