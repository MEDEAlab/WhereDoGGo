#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script will download contigs (each genome's latest version) from Genbank for a list of assemblies. It is essentially a wrapper for ncbi-datasets-cli with small extra features (e.g., an extra failed/non-genome download check).

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#1) ncbi-datasets-cli (https://github.com/ncbi/datasets)

assemblies="$1"

#Check if all dependencies are installed.
if ! command -v datasets &> /dev/null
then
    echo "Program datasets not installed. Download it from https://github.com/ncbi/datasets. Exiting."
    exit 1
fi

cat << EndOfMessage
#Script: downloadcontigs.sh
#Version: v20240713
#Usage: downloadcontigs.sh <assemblies>
#<assemblies> must be a text file of versionless Genbank assemblies (1/line). RefSeq will be converted to Genbank and version numbers removed. (required)
#For more information refer to the comments in the script and/or the Github page.
EndOfMessage

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
#To make it path and extension agnostic (since it's input and not created by another command to be certain of the name), basename removes the path and the perl one-liner the extensions.
#Output is always in the working directory.
baseassemblies="$(basename "$assemblies" | perl -p -e 's/^(.*?)\..*/$1/g')"

#For any given assemblies file, remove all the files a previous run would have produced to avoid clashes.
#All output will be in the working directory.
echo "Removing files and directories with names identical to the output."
rm -r "$baseassemblies"_contigs/ "$baseassemblies".failed "$baseassemblies".retry ncbi_dataset.zip ncbi_dataset/ 2> /dev/null

creationtime="$(echo $(TZ=UTC date '+%Y-%m-%d-%H-%M-%S'))"
echo "This script was run at $creationtime UTC."

#Fix the formatting of the assemblies file. RefSeq gets converted to Genbank and version numbers are removed. Within WhereDoGGo? the formatting should already be correct.
#The format for Genbank assemblies is GCA_NNNNNNNNN.N where N is [0-9]. The final .N is the version number. Refseq assemblies start with "GCF" instead.
#There's an assumption that the assemblies file already contains assembly accessions (at least) and not sth entirely different.
echo "Converting assemblies to versionless Genbank."
#Convert RefSeq assemblies to Genbank.
perl -p -i -e 's/^GCF_/GCA_/g' "$assemblies"
#Remove the assembly version.
perl -p -i -e 's/^(GCA_\d+?)\.\d+/$1/g' "$assemblies"

#Download genomes (Genbank, latest version) with ncbi-datasets-cli.
echo "Downloading contig datasets from GenBank."
datasets download genome accession --inputfile "$assemblies" --assembly-source GenBank --assembly-version latest

#Unzip the contig fna files from the downloaded zip file (we haven't changed the default name ncbi_dataset.zip). We don't care about the readme, json assembly reports etc. Then remove the zip file to save disk space.
echo "Uncompressing genome files."
touch "$baseassemblies".retry
mkdir "$baseassemblies"_contigs/
#Unzip the fna files.
#If any of the downloaded genome files are corrupted (usually because of an unstable internet connection), an error message "error:  invalid compressed data to inflate" followed by the file will be passed.
#From these messages we extract versionless assemblies and redirect them to the .retry file, ignoring other error messages.
unzip -q -j -d "$baseassemblies"_contigs/ ncbi_dataset.zip '*.fna' 2>&1 >/dev/null | grep 'invalid compressed data to inflate' | perl -p -e 's/^.*?_contigs\/(GCA_\d+?)\..*/$1/g' > "$baseassemblies".retry
rm ncbi_dataset.zip

#Create a variable that will be used as a checkpoint for any downloaded corrupted genomes.
declare -i corruptedcheck
corruptedcheck=0

#If the .retry file is not empty, we enter a loop of retrying the downloads of any corrupted genomes until there are none left.
if [ "$(wc -l < "$baseassemblies".retry)" -gt 0 ]
then
	#Check the value of the corruptedcheck variable.
	while [ "$corruptedcheck" -eq 0 ] ; do
		#Check how many lines the .retry file contains, to pass a grammatically correct message.
		if [ "$(wc -l < "$baseassemblies".retry)" -eq 1 ]
		then
			echo "One corrupted genome file detected in the previous download. Retrying for that genome."
		else
			echo ""$(wc -l < "$baseassemblies".retry)" corrupted genome files detected in the previous download. Retrying for these genomes."
		fi
		#Remove the corrupted files to avoid conflicts with the new download.
		echo "Removing corrupted genome files."
		while IFS= read -r line; do rm -r "$baseassemblies"_contigs/"$line"* ; done < "$baseassemblies".retry
		#Reattempt the downloads using the assembly list in the .retry file.
		echo "Downloading contig datasets from GenBank."
		datasets download genome accession --inputfile "$baseassemblies".retry --assembly-source GenBank --assembly-version latest
		echo "Uncompressing genome files."
		unzip -q -j -d "$baseassemblies"_contigs/ ncbi_dataset.zip '*.fna' 2>&1 >/dev/null | grep 'invalid compressed data to inflate' | perl -p -e 's/^.*?_contigs\/(GCA_\d+?)\..*/$1/g' > "$baseassemblies".retry
		rm ncbi_dataset.zip
		#If the retry file is now empty, change the value of corruptedcheck to break the loop.
		if [ "$(wc -l < "$baseassemblies".retry)" -eq 0 ]
		then
			corruptedcheck=1
			echo "No corrupted genome files detected in the previous download. Proceeding."
			rm "$baseassemblies".retry
		fi
	done
else
	echo "No corrupted genome files detected in the previous download. Proceeding."
	rm "$baseassemblies".retry
fi

#Rename the dataset name stems to only a versionless assembly. This will make it easier to modify the headers later.
#The new name is the minimum amount of characters until we encounter the first dot, including the preceding directory structure. We substitute the remaining characters with the extension .fna .
echo "Modifying genome filename stems to correspond to the versionless assemblies."
#This renaming command works too but it seems a bit sketchier.
#for i in "$baseassemblies"_contigs/*.fna ; do mv "$i" "${i%%.*}".fna ; done
for i in "$baseassemblies"_contigs/*.fna ; do mv "$i" "$baseassemblies"_contigs/"$(basename "$i" | perl -p -e 's/^(.*?)\..*/$1/g')".fna ; done

#Check if any of the downloads failed or shouldn't be considered a genome, as it contains very little data. The genome afiles are removed and their assemblies logged in the failed file.
#We've arbitrarily selected a size of 10 bytes as the threshold.
echo "Removing failed or non-genome downloads and writing assemblies to the .failed file."
touch "$baseassemblies".failed
for i in "$baseassemblies"_contigs/*.fna ; do
	if [ $(wc -c < $i) -lt 10 ] ; then
		echo "$(basename "$i" .fna)" >> "$baseassemblies".failed
		rm "$i"
	fi
done

#Check if all the asssemblies were downloaded and are still there. If not, exit with a warning code. doggo_fetch will exit if it sees that code.
declare -i dlwarn
dlwarn=0
while IFS= read -r line; do
	if [ "$(ls "$baseassemblies"_contigs/"$line".fna | wc -l)" -eq 0 ] ; then
    echo "$line" >> "$baseassemblies".failed
		dlwarn=5
	fi
done < "$assemblies"

if [ "$dlwarn" -eq 5 ]
then
	echo "All done! WARNING: ONE OR MORE GENOMES ARE MISSING (CHECK THE LOG AND FAILED FILES)."
	exit 5
fi

#Congrats, you're done!
echo "All done!"
