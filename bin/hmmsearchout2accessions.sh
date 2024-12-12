#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script extracts dereplicated sequence accessions from hmmsearch's --tblout or --domtblout output file.

#Dependencies
#NONE

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubunti. Please report any issues.
#NOTE 2: Accessions are dereplicated because hmmsearch can get multiple hits in the same sequence for a given domain in --domtblout. For --tblout there's no difference either way. TODO: Test how seqtk handles multiple instances of the same accession.

hmmsearchouts="$1"
filext="$2"

cat << EndOfMessage
#Script: hmmsearchout2accessions.sh
#Version: v20241212
#Usage: hmmsearchout2accessions.sh <hmmsearchouts> <filext>
#<hmmsearchouts> must be the path to the directory containing the hmmsearch output files. (trailing slash optional) (required)
#<filext> must be the filename extension of the HMMER -tblout files from which accessions will be extracted. (leading dot optional) (required)
#Accessions must not contain spaces, since it's used as a delimiter in the file.
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

#Check if the directory containing the hmmsearch output files exists, otherwise exit.
if [ -d $hmmsearchouts ]
then
	echo "Contig directory found. Proceeding."
else
	echo "Contig directory not found. Exiting."
	exit 1
fi

#Check if the user provided a path ending with a slash ("/"), otherwise add it.
if [[ $hmmsearchouts != *\/ ]]
then
	hmmsearchouts="${hmmsearchouts}/"
fi

#Check if the extension provided starts with a dot, otherwise add it.
if [[ $filext != \.* ]]
then
	filext=".${filext}"
fi

#Check if there exists at least one file with the chosen extension in the working directory, otherwise exit.
#As opposed to the version of this find command that includes $contigs, we can't put quotes around the dot or declare it as a variable name.
if [ $(find "$hmmsearchouts" -mindepth 1 -maxdepth 1 -name "*$filext" | wc -l | sed 's/ //g') -gt 0 ]; then
  echo "File(s) with the given extension found in the hmmsearchouts directory. Proceeding."
else
  echo "No files with given extension found in the hmmsearchouts directory. Exiting."
  exit 1
fi

#Remove files from previous runs.
echo "Removing files with names identical to the output."
rm -r *.accessions 2> /dev/null

#Extract accessions from HMMER outfile.
echo "Extracting accessions."
for i in "${hmmsearchouts}"*"${filext}" ; do perl -n -e 'print if not m/^#.*/g' $i | perl -p -e 's/^(.*?) .*/$1/g' | sort | uniq >> "$(basename "$i" "${filext}")".accessions ; done

#Check if any searches resulted in no hits, that is empty .accessions files.
echo "Removing markers with no hits."
for i in *.accessions ; do
	if [ "$(wc -l < "$i" | sed 's/ //g')" -eq 0 ] ; then
		echo "WARNING: "$i" is an empty file. Probably no hits found by hmmsearch. Removing."
		rm -r "$i"
	elif [ "$(wc -l < "$i" | sed 's/ //g')" -eq 1 ] ; then
		echo "WARNING: "$i" contains one accession. Aligning and trimming is impossible (and it might bias the phylogenies). Removing."
		rm -r "$i"
	fi
done

#Congrats, you're done!
echo "All done!"
