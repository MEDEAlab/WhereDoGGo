#!/bin/bash

##The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script outputs the taxa names, dataset filenames, number of taxa per dataset, number of markers/sequences per taxon, and number of markers/sequences per taxon per dataset (also checking for multiplicates) for all files with a given extension in the working directory.

#Dependencies
#NONE

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.

alignments="$1"
filext="$2"

cat << EndOfMessage
#Script: preconcatenation.sh
#Version: v20241212
#Usage: preconcatenation.sh <alignments> <filext>
#<alignments> must be the path to the directory containing the alignment files. (trailing slash optional) (required)
#<filext> must be the filename extension of the alignment files. (leading dot optional) (required)
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

#Check if the directory containing the alignments exists, otherwise exit.
if [ -d $alignments ]
then
	echo "Alignments directory found. Proceeding."
else
	echo "Alignments directory not found. Exiting."
	exit 1
fi

#Check if the user provided a path ending with a slash ("/"), otherwise add it.
if [[ $alignments != *\/ ]]
then
	alignments="${alignments}/"
fi

#Check if the extension provided starts with a dot, otherwise add it.
if [[ $filext != \.* ]]
then
	filext=".${filext}"
fi

#Check if there exists at least one file with the chosen extension in the working directory, otherwise exit.
if [ $(find "$alignments" -mindepth 1 -maxdepth 1 -name "*$filext" | wc -l | sed 's/ //g') -gt 0 ]; then
  echo "File(s) with the given extension found in the alignments directory. Proceeding."
else
  echo "No files with given extension found in the alignments directory. Exiting."
  exit 1
fi

#Remove files from previous runs.
echo "Removing files with names identical to the output."
rm -r taxa.names dataset.names dataset.numbers taxa.numbers multiples.check dataset.pass dataset.notpass 2> /dev/null

#Determine the number of taxa that will go into the concatenation.
echo "Creating taxa.names file."
for i in "${alignments}"*"${filext}" ; do grep ">" $i | perl -p -e 's/>(.*?) .*/$1/gm' ; done | sort | uniq >> taxa.names

#Write the names of the datasets to be concatenated.
echo "Creating dataset.names file."
for i in "${alignments}"*"${filext}" ; do echo "$(basename $i)" >> dataset.names ; done

#Calculate the number of taxa per file.
echo "Creating dataset.numbers file."
for i in "${alignments}"*"${filext}" ; do
  count1=$(grep ">" "$i" | wc -l | sed 's/ //g')
	dataset1=$(basename $i)
  echo -e "${dataset1}\t${count1}" >> dataset.numbers
done

#Calculate the number of sequences per taxon.
echo "Creating taxa.numbers file."
while IFS= read -r line ; do echo -e "$line\t$(egrep "^>$line " "${alignments}"*"${filext}" 2> /dev/null | wc -l | sed 's/ //g')" >> taxa.numbers ; done < taxa.names

#Check for multiple sequences.
echo "Checking for multiplicate sequences. Creating multiples.check file."
while IFS= read -r line; do
  for i in "${alignments}"*"${filext}"; do
    count2=$(egrep "^>$line " "$i" 2> /dev/null | wc -l | sed 's/ //g')
		dataset2=$(basename $i)
    echo -e "${dataset2}\t${line}\t${count2}" >> multiples.check
  done
done < taxa.names

while IFS= read -r line; do
	declare -i multicheck
	multicheck=$(echo "$line" | perl -p -e 's/^.*?\t.*?\t(.*?)/$1/g')
	if [ "$multicheck" -gt 1 ]
	then
		echo "Taxa with more than one sequence in a dataset found. Exiting."
		exit 1
	fi
done < multiples.check


#Check which datasets have >50% of the taxa; the rest will be ignored for the final concatenation.
echo "Checking which datasets have >50% of taxa present. Creating dataset.pass and dataset.notpass files."
touch dataset.pass dataset.notpass
threshold="$(($(wc -l < taxa.names | sed 's/ //g')/2))"
while IFS= read -r line; do
  threshcheck=$(echo "$line" | perl -p -e 's/^.*?\t(\d+)/$1/')
  whichdataset=$(echo "$line" | perl -p -e 's/^(.*?)\t.*/$1/')
  if python -c "import sys; print(float(sys.argv[1]) > float(sys.argv[2]))" "$threshcheck" "$threshold" | grep -q True
	then
    echo "$whichdataset" >> dataset.pass
  else
    echo "$whichdataset" >> dataset.notpass
  fi
done < dataset.numbers

#Congrats, you're done!
echo "All done!"
