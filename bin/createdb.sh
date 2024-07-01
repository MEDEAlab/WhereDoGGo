#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script will cat orfs produced by contigs2orfs into a single database against which we can then run homology searches.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.
#NOTE 2: This is the WhereDoGGo? version of the script that assumes that a properly named directory containing the orf faa files is present in the working directory (normally created by contigs2orfs.sh).
#All faa files in the orfs directory are included in the db indiscriminately.
#TODO: Have the script run off a list of assemblies instead (to use subsets of the genomes)?

#Dependencies
#NONE

assemblies="$1"

cat << EndOfMessage
#Script: createdb.sh
#Version: v20240701
#Usage: createdb.sh <assemblies>
#<assemblies> must be a text file of versionless assemblies (1/line). (required)
#For more information refer to the comments in the script and/or the Github page.
EndOfMessage
#We don't really care about the contents of the assemblies file here, just its filename stem for finding the orf directory. Thus, there are no checks for proper formatting of the assembly list.
#TODO: Include an option for picking which directory contains the orf files (and maybe even their file extension).

#Check if the number of arguments is correct, otherwise exit.
if [ "$#" -eq 1 ]
then
	echo "One argument found. Proceeding."
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
if [ -d "$baseassemblies"_orfs ]
then
	echo "ORF directory found. Proceeding."
else
	echo "ORF directory not found in the working directory. Exiting."
	exit 1
fi

#Remove the output from any previous runs to avoid appending.
echo "Removing files with names identical to the output."
rm -r "$baseassemblies".database 2> /dev/null

#Cat all the orf files together.
echo "Creating local database."
for i in "$baseassemblies"_orfs/*.faa ; do cat "$i" >> "$baseassemblies".database ; done

#Congrats, you're done!
echo "All done!"
