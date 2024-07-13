#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script will create assemblies (list of assemblies) and assmebliesnames (assemblies and species names, tab-delimited) files for a local set of genome files (as contigs).
#Since these genomes did not originate from NCBI, we have to create dummy assemblies and species names.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#NONE

contigs="$1"
filext="$2"
code="$3"

cat << EndOfMessage
#Script: createrecords.sh
#Version: v20240713
#Usage: createrecords.sh <contigs> <filext> <code>
#<contigs> must be the path to the directory containing the genome files (as contigs). (trailing slash optional) (required)
#<filext> must be the file extension of the genome files. (leading dot optional) (required)
#<code> must be a project code (alphanumeric and/or underscores only) that will be used as a filename stem and in assemblies. If one is not given, it will be randomly generated. (optional)
#For more information refer to the comments in the script and/or the Github page.
EndOfMessage

#Check if the number of arguments is correct, otherwise exit.
if [[ "$#" -eq 2 || "$#" -eq 3 ]]
then
	echo "$# arguments found. Proceeding."
else
	echo "Wrong number of arguments given. Exiting."
	exit 1
fi

#Check if the directory containing the genomes exists, otherwise exit.
if [ -d $contigs ]
then
	echo "Contig directory found. Proceeding."
else
	echo "Contig directory not found. Exiting."
	exit 1
fi

#Check if the user provided a path ending with a slash ("/"), otherwise add it. This relates to finding files later.
if [[ $contigs != *\/ ]]
then
	contigs="${contigs}/"
fi

#Check if the extension provided starts with a dot, otherwise add it. This relates to finding files later.
if [[ $filext != \.* ]]
then
	filext=".${filext}"
fi

#Check if there exists at least one file with the chosen extension in the genome directory, otherwise exit.
if [ $(find "$contigs" -mindepth 1 -maxdepth 1 -name "*$filext" | wc -l) -gt 0 ]; then
  echo "File(s) with the given extension found in the contig directory. Proceeding."
else
  echo "No files with given extension found in the contig directory. Exiting."
  exit 1
fi

#If project code is given, check if it is alphanumeric and/or underscores. Otherwise, randomly generate one consisting of five capital letters.
if [ "$#" -eq 3 ]
then
	if [[ "$code" =~ ^[a-zA-Z0-9_]+$ ]]
	then
		echo "Project code is valid. Proceeding."
	else
		echo "Project code is invalid. Exiting."
		exit 1
	fi
elif [ "$#" -eq 2 ]
then
	code="$(LC_ALL=C tr -dc A-Z </dev/urandom | head -c 5)"
	echo "Generated random project code ${code}. Proceeding."
fi

#Remove the output from any previous runs to avoid appending.
echo "Removing files with names identical to the output."
rm -r "$code".assemblies "$code".assembliesnames 2> /dev/null

#Create a six digit counter for the genome dummy assemblies.
counter=100001

#Create the assemblies files, where each assembly is a fusion of the project code with the iterating counter.
#Also create the tab-delimited assembliesnames file where the species name is the stem of each genome filename.
echo "Creating the .assemblies and .assembliesnames files."
for i in "${contigs}"*"${filext}" ; do
	echo -e "${code}${counter}\t$(basename "$i" "${filext}")" >> "$code".assembliesnames
	echo  "${code}${counter}" >> "$code".assemblies
	counter=$((counter+1))
done

#Congrats, you're done!
echo "All done!"
