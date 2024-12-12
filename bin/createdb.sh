#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script will cat orfs produced by contigs2orfs into a single database against which we can then run homology searches.

#NOTE 1: All code was written and tested on or ARM Intel macOS and Ubuntu. Please report any issues.
#All files in the orfs directory with a given extension are included in the db indiscriminately.
#TODO: Have the script run off a list of assemblies instead (to use subsets of the genomes)?

#Dependencies
#NONE

assemblies="$1"
orfs="$2"
filext="$3"

cat << EndOfMessage
#Script: createdb.sh
#Version: v20241212
#Usage: createdb.sh <assemblies> <orfs> <filext>
#<assemblies> must be a text file of versionless assemblies (1/line). (required)
#<orfs> must be the path to the directory containing the genome files (as orfs). (trailing slash optional) (required)
#<filext> must be the extension of the genome files. (leading dot optional) (required)
#For more information refer to the comments in the script and/or the Github page.
EndOfMessage
#We don't really care about the contents of the assemblies file here, just its filename stem for finding the orf directory. Thus, there are no checks for proper formatting of the assembly list.

#Check if the number of arguments is correct, otherwise exit.
if [ "$#" -eq 3 ]
then
	echo "Three arguments found. Proceeding."
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

#Isolate the stem of the assemblies file name as a separate variable to avoid any issues with extensions.
baseassemblies="$(basename "$assemblies" | perl -p -e 's/^(.*?)\..*/$1/g')"

#Check if the orfs directory exists, otherwise exit.
if [ -d $orfs ]
then
	echo "ORF directory found. Proceeding."
else
	echo "ORF directory not found. Exiting."
	exit 1
fi

#Check if the user provided a path ending with a slash ("/"), otherwise add it. This relates to finding files later.
if [[ $orfs != *\/ ]]
then
	orfs="${orfs}/"
fi

#Check if the extension provided starts with a dot, otherwise add it. This relates to finding files later.
if [[ $filext != \.* ]]
then
	filext=".${filext}"
fi

#Check if there exists at least one file with the chosen extension in the orfs directory, otherwise exit.
#If this find command is run without $orfs assigned, it will fail (without the quotes it will default to current directory). If $filext is not assigned it will find everything (because it defaults to an asterisk).
if [ $(find "$orfs" -mindepth 1 -maxdepth 1 -name "*$filext" | wc -l | sed 's/ //g') -gt 0 ]; then
  echo "File(s) with the given extension found in the ORF directory. Proceeding."
else
  echo "No files with given extension found in the ORF directory. Exiting."
  exit 1
fi

#Remove the output from any previous runs to avoid appending.
echo "Removing files with names identical to the output."
rm -r "$baseassemblies".database 2> /dev/null

#Cat all the orf files together.
echo "Creating local database."
for i in "${orfs}"*"${filext}" ; do cat "$i" >> "$baseassemblies".database ; done

#Congrats, you're done!
echo "All done!"
