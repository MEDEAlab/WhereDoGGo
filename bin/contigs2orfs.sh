#!/bin/bash

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script will predict ORFs with Pyrodigal for the contig/scaffold genome files in a directory (names given as a list of assemblies, file extension of the files is given by the user)
#TODO: Make another version that uses Bakta (or Prokka) to include annotations in the headers.

#NOTE 1: All code was written and tested on Intel or ARM macOS and Ubuntu. Please report any issues.
#NOTE 2: This is the WhereDoGGo? version of the script that assumes that a correctly formatted (genbank only, versionless, 1/line, assembly & name tab-delimited) assembliesnames file is present in the working directory.
#All genome files in the contigs directory are used for ORF prediction indiscriminately.
#TODO: Have the script run off a list of assemblies instead?

#Dependencies
#1) Pyrodigal (https://github.com/althonos/pyrodigal)

assemblies="$1"
contigs="$2"
filext="$3"
whatwematch="$4"
code25="$5"

#Check if all dependencies are installed.
if ! command -v pyrodigal &> /dev/null
then
    echo "Program pyrodigal not installed. Download it from https://github.com/althonos/pyrodigal. Exiting."
    exit 1
fi

cat << EndOfMessage
#Script: contigs2orfs.sh
#Version: v20241212
#Usage: contigs2orfs.sh <assemblies> <contigs> <filext> <whatwematch> <code25>
#<assemblies> must be a text file of versionless assemblies (1/line). (required)
#<contigs> must be the path to the directory containing the genome files (as contigs). (trailing slash optional) (required)
#<filext> must be the extension of the genome files. (leading dot optional) (required)
#<whatwematch> must be "assembies" or "names". (case-sensitive) (required)
#<code25> must be "code25". If added, the pyrodigal runs will default to using Genetic Code 25. (optional)
#For more information refer to the comments in the script and/or the Github page.
EndOfMessage
#We don't really care about the contents of the assemblies file here, just its filename stem for naming. Thus, there are no checks for proper formatting of the assembly list.

#Check if the number of arguments is correct, otherwise exit.
if [[ "$#" -eq 4 || "$#" -eq 5 ]]
then
	echo "$# arguments found. Proceeding."
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

#Check if the .assembliesnames file exists in the working directory, otherwise exit.
#TODO: Maybe check if the assembliesnames content and formatting is correct? (and thus remove NOTE 2)
if [ -f "$baseassemblies".assembliesnames ]
then
	echo "Assembliesnames file found. Proceeding."
else
	echo "Assembliesnames file not found in the working directory. Exiting."
	exit 1
fi

#Check if the contigs directory exists, otherwise exit.
if [ -d $contigs ]
then
	echo "Contig directory found. Proceeding."
else
	echo "Contig directory not found. Exiting."
	exit 1
fi

#Check if the whatwematch variable is correct (case-sensitive).
if [[ "$whatwematch" = "assemblies" || "$whatwematch" = "names" ]] ; then
	echo "Matching "$whatwematch" in assembliesnames file. Proceeding."
else
	echo "Invalid whatwematch option. Exiting."
	exit 1
fi

#Check if code25 argument is code25
if [ "$#" -eq 5 ]
then
  if [ "$code25" = "code25" ]
  then
    echo "Code25 argument found. Proceeding."
  else
    echo "Fifth argument is not code25. Exiting."
    exit 1
  fi
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
#We were checking with ls before but when trying to download all Bacteria, it goes way above ARGMAX.
#If this find command is run without $contigs assigned, it will fail (without the quotes it will default to current directory). If $filext is not assigned it will find everything (because it defaults to an asterisk).
if [ $(find "$contigs" -mindepth 1 -maxdepth 1 -name "*$filext" | wc -l | sed 's/ //g') -gt 0 ]; then
  echo "File(s) with the given extension found in the contig directory. Proceeding."
else
  echo "No files with given extension found in the contig directory. Exiting."
  exit 1
fi

#For any given assemblies file, remove all the files a previous run would have produced to avoid clashes.
echo "Removing files and directories with names identical to the output."
rm -r "$baseassemblies"_orfs/ "$baseassemblies".pyrodigallog 2> /dev/null

#Predict ORFs with Pyrodigal.
echo "Predicting ORFs with Pyrodigal."
#Create the orfs directory. Pyrodigal requires all directory structures to already be in place.
mkdir "$baseassemblies"_orfs
#Create a directory to retain the pyrodigal output files.
mkdir "$baseassemblies"_orfs/pyrodigal_runs
#We run Pyrodigal on all genomes in the contigs directory indiscriminately.
for i in "${contigs}"*"${filext}" ; do
	#Create the assemblyacc variable for filenames and logging. This is the assembly or name corresponding to each genome file, as determined by the whatwematch variable.
	filestem="$(basename "$i" "${filext}")"
	if [ "$whatwematch" = "assemblies" ] ; then
		assemblyacc="$(perl -n -e 'BEGIN { $matcher = pop } print if m/^$matcher\t.*/' "$baseassemblies".assembliesnames "$filestem" | perl -p -e 's/^(.*?)\t.*/$1/g')"
	elif [ "$whatwematch" = "names" ] ; then
		assemblyacc="$(perl -n -e 'BEGIN { $matcher = pop } print if m/^.*?\t$matcher\n/' "$baseassemblies".assembliesnames "$filestem" | perl -p -e 's/^.*?\t(.*)/$1/g')"
	fi
	#Pyrodigal requires all directory structures to be already in place, so create individual directories for the output of each genome.
	mkdir "$baseassemblies"_orfs/"$assemblyacc"
	#Echo the assembly into the log to be able to distinguish outputs.
	echo "$assemblyacc" >> "$baseassemblies".pyrodigallog
  #echo "$assemblyacc"
	#Run Pyrodigal on each contig file, and output its faa (aa ORFs), gbk (old Genbank format), and ffn (nucleotide ORFs) files. Redirect both output and errors to the log.
	#At NCBI, gbk has been substituted by gbff but it's not an option, plus some programs take gbk as input. Alternatively, change the output to gff.
	#Unlike Prodigal, which iirc freaks out if you -o to a gbk file, since it does this by default and will only take gff output, Pyrodigal needs this argument otherwise it outputs the gbk contents to stdout (or pyrodigallog here, since we redirect).
	#For the Pyrodigal run, if the genome is from SR1 or Gracilibacteria (c__JAEDAM01), use Code 25.
	if grep "$filestem" "$baseassemblies".assembliesnames | grep -q "c__JAEDAM01"
	then
		(pyrodigal -m -i "$i" -a "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".faa -o "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".gbk -d "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".ffn -g 25 ) >> "$baseassemblies".pyrodigallog 2>&1
	elif [ "$#" -eq 5 ]
	then
		(pyrodigal -m -i "$i" -a "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".faa -o "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".gbk -d "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".ffn -g 25 ) >> "$baseassemblies".pyrodigallog 2>&1
	else
		(pyrodigal -m -i "$i" -a "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".faa -o "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".gbk -d "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".ffn ) >> "$baseassemblies".pyrodigallog 2>&1
	fi
	#Echo a line with two slashes in the log after each output. This is a commmon entry delimiter e.g., for Uniprot files.
	echo "//"  >> "$baseassemblies".pyrodigallog
	#Create a copy of the faa file to be used for header modification and asterisk removal.
	cp "$baseassemblies"_orfs/"$assemblyacc"/"$assemblyacc".faa "$baseassemblies"_orfs/
	#Move the pyrodigal run output with all the runs.
	mv "$baseassemblies"_orfs/"$assemblyacc" "$baseassemblies"_orfs/pyrodigal_runs/
	done

#Edit the fasta headers.
echo "Editing faa headers and removing asterisks corresponding to stop codons at the end of sequences."
for i in "$baseassemblies"_orfs/*.faa ; do
	#Redefine the filestem variable for faa files (that now corresponds to the assembly, like assemblyacc did previously). This makes it easier to define the modder variable.
	filestem="$(basename "$i" .faa)"
	#Define the modder variable, the line from assembliesnames containing the genome's assembly (and species name), depending on the whatwematch variable (we're trying to prevent substring matches).
	if [ "$whatwematch" = "assemblies" ] ; then
		modder="$(perl -n -e 'BEGIN { $matcher = pop } print if m/^$matcher\t.*/' "$baseassemblies".assembliesnames "$filestem")"
	elif [ "$whatwematch" = "names" ] ; then
		modder="$(perl -n -e 'BEGIN { $matcher = pop } print if m/^.*?\t$matcher\n/' "$baseassemblies".assembliesnames "$filestem")"
	fi
	#Remove all asterisks from the faa files. These (as far as we've encountered) correspond to stop codons.
	#TODO: At some point we need to check how often stop codons can occur mid-sequence by checking all Bacteria and Archaea. We don't really care about them in headers (and in WhereDoGGo? they shouldn't occur anyway).
	perl -p -i -e 's/\*//g' "$i"
	#Modify the pyrodigal output faa headers by keeping only the first element (accession kinda), and adding the modder separated by a tab.
	#We get rid of the rest of the header as it's not particularly informative for WhereDoGGo?.
	perl -p -i -e 'BEGIN { $matcher = pop } s/^>(.*?) .*/>$1\t$matcher/g' "$i" "$modder"
	#Rearrange the fasta headers to a pseudo-NCBI format aka the doggo format (>accession assembly [species]). The assembly is tagged on as part of the accesions to ensure their uniqueness.
	#TODO: Reconsider these header edits? Perhaps don't get rid of the pyrodigal orf info entirely.
	perl -p -i -e 's/>(.*?)\t(.*?)\t(.*)/>$2_$1 $2 \[$3\]/g' "$i"
done

#Congrats, you're done!
echo "All done!"
