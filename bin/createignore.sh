#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script checks which assemblies in from GTDB metadata are not found among all NCBI assemblies for a given domain. These become the ignore list. They are also looked up against NCBI for their organism names.
#The process is repeated for atypical assemblies in NCBI and local asssemblies neither in the non-atypical list, nor in the ignore list. For those their names and atypical warnings are also fetched.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.
#TODO: In its current version, the script misses some edge cases e.g., ones that have been removed/suppressed from the Assembly database but a genome is available in the Nucleotide database. Should update with double-checking through full domain downloads.

#Dependencies
#1) ncbi-datasets-cli (https://github.com/ncbi/datasets)

input_file="$1"
taxon="$2"
domain="$3"

#Check if all dependencies are installed.
if ! command -v datasets &> /dev/null
then
    echo "Program datasets not installed. Download it from https://github.com/ncbi/datasets. Exiting."
    exit 1
fi

cat << EndOfMessage
#Script: createignore.sh
#Version: v20240701
#Usage: createignore.sh <input_file> <taxon> <domain>
#<input_file> must be a file with tab-delimited GTDB metadata (original or parsed). (required)
#<taxon> must be the GTDB taxon for which the ignorelist will be created. (required)
#<domain> must be "Bacteria" or "Archaea". (case-sensitive) (required)
#For more information refer to the comments in the script and/or the Github page.
EndOfMessage

#Check if the number of arguments is correct, otherwise exit.
if [ "$#" -eq 3 ]
then
	echo "Three arguments found. Proceeding."
else
	echo "Wrong number of arguments given. Exiting."
	exit 1
fi

#Check if the assemblies file exists, otherwise exit.
if [ -f $input_file ]
then
	echo "Input metadata file found. Proceeding."
else
	echo "Input metadata file not found. Exiting."
	exit 1
fi

#Check if the domain name is correct (case-sensitive).
if [ "$domain" = "Bacteria" ] ; then
	echo "Domain Bacteria found. Proceeding."
	domaintxid="2"
elif [ "$domain" = "Archaea" ] ; then
	echo "Domain Archaea found. Proceeding."
	domaintxid="2157"
else
	echo "Invalid domain name. Exiting."
	exit 1
fi

#Create the list of assemblies. If it's empty (taxon doesn't exist), remove it and exit.
#The conversion to versionless Genbank is included.
echo "Removing any existing assemblies files for the taxon given and creating a new one."
rm -r "$taxon".assemblies 2> /dev/null
perl -p -e 's/^\w*?_(.*?)\.\d+\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t.*?\t(.*?)\t.*/$1\t$2/g' "$input_file" | grep "$taxon" | perl -p -e 's/^(.*?)\t.*/$1/g' | perl -p -e 's/^GCF_/GCA_/g' >> "$taxon".assemblies
if [ $(wc -l < "$taxon".assemblies) -gt 0 ] ; then
  echo "Assemblies found for the taxon given in the input metadata file. Proceeding."
else
  echo "No assemblies found for the taxon given in the input metadata file (resulted in empty file). Exiting."
  rm -r "$taxon".assemblies
  exit 1
fi

#For any given assemblies file, remove all the files a previous run would have produced to avoid clashes.
#All output will be in the working directory.
echo "Removing all existing output files for the same domain."
rm -r "$domain"*_NCBI.all "$domain"*.ignore "$domain"*.ignorenames "$domain"*.ignorediff "$domain"*_NCBI.noatypical "$domain"*.atypical "$domain"*.atypicalnameswarnings "$domain"*.atypicaldiff 2> /dev/null

creationtime="$(echo $(TZ=UTC date '+%Y-%m-%d-%H-%M-%S'))"
echo "This script was run at $creationtime UTC."

#Create a file with all the assemblies in NCBI for the selected domain. Suppressed assemblies are excluded automatically (but it probably doesn't matter, since we go for Genbank).
#TODO: Previously the command also included --exlude-atypical but this removes too many species for Archaea. Maybe add it back to be on the safe side?
echo "Creating a file with all assemblies in NCBI for $domain."
datasets summary genome taxon $domaintxid --as-json-lines --assembly-source genbank --assembly-version latest | dataformat tsv genome --fields accession --elide-header | perl -p -e 's/^(.*?)\.\d+/$1/g' >> "$domain"_"$creationtime"_NCBI.all

#Check which of the assemblies in the assemblies files are found among the NCBI assemblies for the domain. Anything not found is added to the ignorelist.
echo "Writing all assemblies not among the NCBI assemblies to the .ignore file."
touch "$domain"_"$creationtime".ignore
awk 'NR==FNR{set[$0]=$0}NR!=FNR{if($0 in set==0){print $0}}' "$domain"_"$creationtime"_NCBI.all "$taxon".assemblies >> "$domain"_"$creationtime".ignore
#while IFS= read -r line ; do if LC_ALL=C fgrep -q -- "$line" "$domain"_"$creationtime"_NCBI.all ; then : ; else echo "$line" >> "$domain"_"$creationtime".ignore ; fi ; done < "$taxon".assemblies

#For all assemblies in the ignorelist, look up their organism name on NCBI and write it together with a versionless assembly. This is to confirm the most common reason for not finding an assembly (especially for Archaea), reassignment to a different domain.
echo "Writing ignore assemblies and the NCBI organism names to the .ignorenames file."
while IFS= read -r line ; do datasets summary genome accession "$line" --as-json-lines --assembly-source GenBank --assembly-version latest | dataformat tsv genome --fields accession,organism-name --elide-header | perl -p -e 's/^(.*?)\.\d+/$1/g' >> "$domain"_"$creationtime".ignorenames ; done < "$domain"_"$creationtime".ignore

#Check for discrepancies between the assemblies in the ignore and the ignorenames files. Output any in the ignorediff file.
echo "Checking for discrepancies between the assemblies in the .ignore and .ignorenames files. If there are any, writing them to the .ignorediff file. These are usually suppressed GenBank assemblies."
touch "$domain"_"$creationtime".ignorediff
#TODO: Fix the awk one-liner.
#awk 'NR==FNR{set[$0]=$0}NR!=FNR{if($0 in set==0){print $0}}' "$domain"_"$creationtime".ignorenames "$domain"_"$creationtime".ignore >> "$domain"_"$creationtime".ignorediff
while IFS= read -r line ; do if LC_ALL=C fgrep -q -- "$line" "$domain"_"$creationtime".ignorenames ; then : ; else echo "$line" >> "$domain"_"$creationtime".ignorediff ; fi ; done < "$domain"_"$creationtime".ignore

#Create a file with all the assemblies in NCBI for the selected domain excluding atypical.
#TODO: It looks like Refseq suppressed (at least some of them) fall under atypical (in Genbank and/or Refseq). Need to confirm to what extent this happens and if the assminfo-suppression-reason field needs to be included.
echo "Creating a file with all assemblies in NCBI for $domain, excluding atypical."
datasets summary genome taxon $domaintxid --as-json-lines --assembly-source genbank --assembly-version latest --exclude-atypical | dataformat tsv genome --fields accession --elide-header | perl -p -e 's/^(.*?)\.\d+/$1/g' >> "$domain"_"$creationtime"_NCBI.noatypical

#Check which of the assemblies in the assemblies file are found neither in the ignorelist nor among the non-atypical NCBI assemblies.
echo "Writing all assemblies not among the non-atypical NCBI assemblies and not in the .ignore file to the .atypical file."
touch "$domain"_"$creationtime".atypical
awk 'NR==FNR{set[$0]=$0}NR!=FNR{if($0 in set==0){print $0}}' "$domain"_"$creationtime"_NCBI.noatypical "$taxon".assemblies >> "$domain"_"$creationtime".atypical
while IFS= read -r line ; do if LC_ALL=C fgrep -q -- "$line" "$domain"_"$creationtime".ignore ; then perl -p -i -e 'BEGIN { $matcher = pop } s/^$matcher\n//g' "$domain"_"$creationtime".atypical "$line" ; fi ; done < "$domain"_"$creationtime".atypical
#while IFS= read -r line ; do if LC_ALL=C fgrep -q -- "$line" "$domain"_"$creationtime"_NCBI.noatypical "$domain"_"$creationtime".ignore ; then : ; else echo "$line" >> "$domain"_"$creationtime".atypical ; fi ; done < "$taxon".assemblies

#For all atypical assemblies, look up their assminfo-atypicalwarnings and write it in a new file after the assembly (tab-delimited). assminfo-suppression-reason doesn't work here, it's probably a Refseq-only field
echo "Writing atypical assemblies, organism names, and atypical warnings to the .atypicalnameswarnings file."
while IFS= read -r line ; do datasets summary genome accession "$line" --as-json-lines --assembly-source GenBank --assembly-version latest | dataformat tsv genome --fields accession,organism-name,assminfo-atypicalwarnings --elide-header | perl -p -e 's/^(.*?)\.\d+/$1/g' >> "$domain"_"$creationtime".atypicalnameswarnings ; done < "$domain"_"$creationtime".atypical

#Check for discrepancies between the assemblies in the ignore and the ignorenames files. Output any in the ignorediff file.
echo "Checking for discrepancies between the assemblies in the .atypical and .atypicalnameswarnings files. If there are any, writing them to the .atypicaldiff file."
touch "$domain"_"$creationtime".atypicaldiff
#TODO: Fix the awk one-liner.
#awk 'NR==FNR{set[$0]=$0}NR!=FNR{if($0 in set==0){print $0}}' "$domain"_"$creationtime".atypicalnameswarnings "$domain"_"$creationtime".atypical >> "$domain"_"$creationtime".atypicaldiff
while IFS= read -r line ; do if LC_ALL=C fgrep -q -- "$line" "$domain"_"$creationtime".atypicalnameswarnings ; then : ; else echo "$line" >> "$domain"_"$creationtime".atypicaldiff ; fi ; done < "$domain"_"$creationtime".atypical

#Congrats, you're done!
echo "All done!"
