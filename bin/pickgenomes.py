#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script picks genomes according to taxonomic level, resolution, and number using GTDB metadata to pick the best quality genome and maximize the taxonomic range of the genomes picked.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#NONE

import math
import os
import sys

print('#Script: pickgenomes.py')
print('#Version: v20240713')
print('#Usage: python pickgenomes.py <input_tsv> <tax_level> <tax_resolution> <number> <ignore_list>')
print('#<input_tsv> must be a tab-delimited file as per GTDB\'s metadata files for Bacteria or Archaea. (required)')
print('#<tax_level> must be the highest taxonomic level for genomes to be selected as presented in GTDB taxonomy strings e.g. p__Asgardarchaeota (i.e., get genomes from within <tax_level>). (required)')
print('#<tax_resolution> must be a taxonomic level lower than <tax_level>. For each <tax_resolution> in <tax_level>, the script picks <number> genomes (or as many as available). It must be p (phylum), c (class), o (order), f (family), g (genus), or all. "all" picks all genomes for <tax_level>. (required)')
print('#<number> must be the number of genomes to be picked per <tax_resolution> or "all". If <tax_resolution> is "all", <number> must also be "all". (required)')
print('#<ignore_list> must be a text file containing genome assembly accessions (once per line, versionless, e.g. GCA_011362025) that should be not included in the pickgenomes process (optional)')
print('#For more information refer to the comments in the script and/or the Github page.')

ignore_list = list()
acc_score = dict() # key = assembly accession, value = score
acc_taxonomy = dict() # key = assembly accession, value = gtdb taxonomy
resolutions = list() # list with unique resolutions (e.g. families) of what the user asks
unique_resolution_accessions = list() # list with the accessions within every unique (specified) resolution
unique_resolution_scores = list() # list with the scores within every unique (specified) resolution
unique_resolutions_minus_one = list() # list with the resolutions-1 within every unique (specified) resolution
unique_resolution_minus_one_accessions = list() # list with accessions within a unique resolution-1
d=0 # variable to help me determine if the user has requested genus as resolution
unique_resolution_minus_one_scores = list() # list with scores within a unique resolution-1
best_scores_from_every_resolution_minus_one = list() # list with the best score out of every resolution-1 of a specified resolution (used for picking one of each resolution-1 if the resolutions-1 are >= sys.argv[4])
already_printed_accessions = list() # list with the already printed accessions, so as to not print them again
empty_res_minus_one = list() # this list helps us check if the res-1 have all their genomes printed, so we can move to the next ones
sub_acc_score = dict()

# checkpoint for number of arguments
if len(sys.argv) == 5 or len(sys.argv) == 6:
    print(str((len(sys.argv)-1)) + ' arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

file_stem = (str(sys.argv[2]) + '_' + str(sys.argv[3]) + '_' + str(sys.argv[4]))

# checkpoint for input file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Input file found. Proceeding.')
else:
    print('Input file not found. Exiting.')
    sys.exit(1)

#Check if input file is parsed GTDB metadata.
with open(sys.argv[1], 'r') as parsedcheck:
    for line in parsedcheck:
        x = line.split('\t')
        if x[18]=='t':
            pass
        else:
            print('Input file does not contain parsed GTDB metadata. Exiting.')
            sys.exit(1)
print('Input file contains parsed GTDB metadata. Proceeding.')

# checkpoint for taxonomic level
with open(sys.argv[1], 'r') as file:
    content = file.read()
    if str(sys.argv[2]+';') in content:
        print('Taxonomic level is valid. Proceeding.')
    else:
        print('Taxonomic level is invalid. Check spelling. Exiting.')
        sys.exit(1)

# checkpoint for taxonomic resolution
tax_res = ['p', 'c', 'o', 'f', 'g']
if sys.argv[3] in tax_res:
    print('Taxonomic resolution is valid. Proceeding.')
elif sys.argv[3] == 'all':
    print('Taxonomic resolution is valid (all). Proceeding.')
else:
    print('Taxonomic resolution is invalid or missing. Check if it is one of: "p", "c", "o", "f", "g", or "all". Exiting.')
    sys.exit(1)

# checkpoint for resolution and number of genomes both being "all"
if sys.argv[3] == 'all' and sys.argv[4] != 'all':
    print('When <tax_resolution> is "all", <number> must also be "all". Exiting.')
    sys.exit(1)
elif sys.argv[3] != 'all' and sys.argv[4] == 'all':
    print('When <number> is "all", <tax_resolution> must also be "all". Exiting.')
    sys.exit(1)

# checkpoint for number of genomes to be picked
flag = True
try:
    int(sys.argv[4])
except:
    flag = False
if flag and int(sys.argv[4])>0:
    print('Number given is valid. Proceeding.')
elif sys.argv[4] == 'all':
    print('Number given is valid (all). Proceeding.')
else:
    print('Number given is invalid or missing. Number must be a natural number or both <tax_resolution> and <number> must be "all". Exiting.')
    sys.exit(1)

if len(sys.argv) == 6: ## this is to make the ignore_list argument optional
    check_ignore = os.path.isfile(sys.argv[5])
    if check_ignore == True:
        print('Ignore list file found. Proceeding.')
        with open(sys.argv[5], 'r') as ignore_text:
            for line in ignore_text:
                x=line.strip()
                ignore_list.append(x)
    else:
        print('Ignore list file not found. Exiting.')
        sys.exit(1)
else:
        print('No ignore list specified. Proceeding.')

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r ' + file_stem + '.assemblies ' + file_stem + '.assembliesnames 2> /dev/null')
os.system(removal)
#TODO: Fix this? It's pointless anyway. Can't do os.WEXITSTATUS, because even if I suppress the output of rm, it exits with an error code.
#if os.WEXITSTATUS(os.system(removal)) == 1:
#    print('Error when removing prior output files. Exiting.')
#    sys.exit(1)

assemblies = open(file_stem+".assemblies", 'a')
assemblies_names = open(file_stem+".assembliesnames", 'a')

print('Picking genomes.')

with open(sys.argv[1], 'r') as taxa_records:
    for line in taxa_records: #1st step: define genome quality scores
        x = line.split('\t')


        if x[0] != 'accession': #Ignore the headers line
            score = 1*float(x[2])-5*float(x[3])+1*(float(x[3])*(float(x[10])/100))+0.5*math.log(int(x[47]),10)+0*math.log(int(x[16]),10) #genome quality calculation (identical to dRep)
            x[0] = x[0][:5] + 'A' + x[0][6:] if x[0][5] == 'F' else x[0]
            if x[0][3:-2] not in ignore_list: ## this is the step where we exclude all the assemblies from the ignore_list
                acc_score.update( {x[0] : score} ) #key = assembly accession, value = score
                acc_taxonomy.update( {x[0] : x[19]} ) #key = assembly accession, value = taxonomy

    if sys.argv[4] == 'all': #1st possibility: the user wants all the genomes of a specific taxonomic lvl
        for key, value in acc_taxonomy.items():
            y = value.split(';')
            if sys.argv[2] in y and key[3:-2] not in ignore_list: #Take all from <tax_level>, other arguments are essentially ignored.
                assemblies.write(key[3:-2]+'\n')

    else:
        if sys.argv[4] == '1': #2nd possibility: the user wants 1 genome per (specified) taxonomic resolution within a specific taxonomic lvl
            for key, value in acc_taxonomy.items(): #For any given entry
                y = value.split(';') #split into individual GTDB taxonomy strings
                if sys.argv[2] in y: #if the <tax_level> exists in the taxonomy strings
                    for item in y: #for each taxonomy string
                        if item[0] == sys.argv[3] and item not in resolutions: #if the taxonomic rank letter corresponds to the first letter of the GTDB taxonomy string and it's not already in resolutions
                            resolutions.append(item) #add that taxonomic rank to unique resolutions (i.e. all possible <tax_resolution> in <tax_level>

            for item in resolutions: #for each entry in unique resolutions
                for key, value in acc_taxonomy.items():
                    y = value.split(';')
                    if sys.argv[2] in y and item in y: #if both the <tax_level> and the resolution are in the taxonomy string
                        unique_resolution_accessions.append(key) #Add the accession of that entry to unique_resolution_accessions, so the list will contain all accessions of each resolution.

                for accession in unique_resolution_accessions: #For each accession in unique_resolution_accessions
                    for key1, value1 in acc_score.items(): #Look in the acc_score dictionary to find that accession
                        if accession == key1:
                            unique_resolution_scores.append(value1) #Append the accession's score to unique_resolution_scores
                for key2, value2 in acc_score.items(): #For each element in acc_score
                    if max(unique_resolution_scores) == value2 and key2 in unique_resolution_accessions: #If the score of the accessions is the highest in unique_resolution_scores,failsafe for multiplicate scores, not to pick random accession with same score, but pick accession that exists in the unique resolution accessions
                        assemblies.write(key2[3:-2] + '\n') #print that accession
                        break #break to go the next resolution, in case of equal scores
                unique_resolution_accessions.clear() #clear the lists
                unique_resolution_scores.clear()

        elif int(sys.argv[4]) > 1: #3rd possibility: the user wants >1 genome per (specified) taxonomic resolution within a specific taxonomic lvl
            for key, value in acc_taxonomy.items(): #create the list of unique resolutions as before
                y = value.split(';')
                if sys.argv[2] in y:
                    for item in y:
                        if item[0] == sys.argv[3] and item not in resolutions:
                            resolutions.append(item)

            for item in resolutions: #for each element in unique resolutions
                #print(item)
                for key, value in acc_taxonomy.items():
                    y = value.split(';')
                    if sys.argv[2] in y and item in y: #if both the <tax_level> and the resolution are in the taxonomy string
                        for item1 in y: # we make a list with all the unique resolutions-1 so as to examine them later
                            if sys.argv[3] == str(y[1])[0] and y[2] not in unique_resolutions_minus_one: #if <tax_resolutions> equals the first letter of y element 1 (aka equals "p" for phylum) and the next lower taxonomic level isn't in the unique resolutions-1, add it.
                                unique_resolutions_minus_one.append(y[2])
                            elif sys.argv[3] == str(y[2])[0] and y[3] not in unique_resolutions_minus_one: #ditto for class
                                unique_resolutions_minus_one.append(y[3])
                            elif sys.argv[3] == str(y[3])[0] and y[4] not in unique_resolutions_minus_one: #ditto for order
                                unique_resolutions_minus_one.append(y[4])
                            elif sys.argv[3] == str(y[4])[0] and y[5] not in unique_resolutions_minus_one: #ditto for family
                                unique_resolutions_minus_one.append(y[5])
                            elif sys.argv[3] == str(y[5])[0]: # this is the special case when the user asks the n best species of the genera within a specified taxonomic lvl
                                d = d+1 #Add 1 to variable d
                                for key, value in acc_taxonomy.items():
                                    y = value.split(';')
                                    if sys.argv[2] in y and item in y: #if both the <tax_level> and the resolution are in the taxonomy string
                                        unique_resolution_accessions.append(key) #Add the accession of that entry to unique_resolution_accessions, so the list will contain all accessions of each resolution.

                                for accession in unique_resolution_accessions: #For each accession in unique_resolution_accessions
                                    for key1, value1 in acc_score.items(): #Look in the acc_score dictionary to find that accession
                                        if accession == key1:
                                            unique_resolution_scores.append(value1) #Append the accession's score to unique_resolution_scores
                                    unique_resolution_scores.sort(reverse = True) # sort them in descending order

                                if len(unique_resolution_scores) <= int(sys.argv[4]): # this is the case of the genomes within the genus being less than or equal to the specified number in sys.argv[4] (<number>)
                                    for item in unique_resolution_accessions: #Take all accessions
                                        assemblies.write(item[3:-2]+'\n')
##                                    for key4, value4 in acc_score.items():
##                                        for item in unique_resolution_scores:
##                                            if item == value4 and key4:
##                                                print(key4)
##                                                break

                                else: # this is the case of the genomes within the resolution being more than the specified number in sys.argv[4]
                                    for item in unique_resolution_scores[:int(sys.argv[4])]: #From the first until <number> scores
                                        for key5, value5 in acc_score.items(): #Take the accession that have that score.
                                            if item == value5 and key5 in unique_resolution_accessions: #failsafe for multiplicate scores, not to pick random accession with same score, but pick accession that exists in the unique resolution accessions
                                                assemblies.write(key5[3:-2]+'\n')
                                                break #break to avoid multiplicates
                                unique_resolution_accessions.clear() #clear the lists
                                unique_resolution_scores.clear()

                if len(unique_resolutions_minus_one) == 1: # this is the case of having only one resolution-1
                    for key, value in acc_taxonomy.items():
                        y = value.split(';')
                        if sys.argv[2] in y and item in y: #if both the <tax_level> and the resolution are in the taxonomy string
                            unique_resolution_accessions.append(key) #Add the accession of that entry to unique_resolution_accessions, so the list will contain all accessions of each resolution.

                    for accession in unique_resolution_accessions: #For each accession in unique_resolution_accessions
                        for key1, value1 in acc_score.items(): #Look in the acc_score dictionary to find that accession
                            if accession == key1:
                                unique_resolution_scores.append(value1) #Append the accession's score to unique_resolution_scores
                        unique_resolution_scores.sort(reverse = True) # sort them in descending order

                    if len(unique_resolution_scores) <= int(sys.argv[4]): # this is the case of the genomes within the resolution being less than or equal to the specified number in sys.argv[4]
                        for access in unique_resolution_accessions:
                            assemblies.write(access[3:-2]+'\n')
##                        for key4, value4 in acc_score.items(): #For each item in acc_score
##                            for item in unique_resolution_scores: #If the unique_resolution_score matches the item's score
##                                if item == value4:
##                                    print(key4) #print the accession
##                                    break #break (go to next) to avoid multiplicates


                    else: # this is the case of the genomes within the resolution being more than the specified number in sys.argv[4]
                        for item in unique_resolution_scores[:int(sys.argv[4])]:
                            for key5, value5 in acc_score.items():
                                if item == value5 and key5 in unique_resolution_accessions:#failsafe for multiplicate scores, not to pick random accession with same score, but pick accession that exists in the unique resolution accessions
                                    assemblies.write(key5[3:-2]+'\n')
                                    break

                    unique_resolution_accessions.clear() #Clear the lists
                    unique_resolution_scores.clear()



                elif d == 0: # which means that the user hasn't asked for species within genera of specified taxonomic lvl AND that the resolutions-1 are more than one within the specified resolution

                    if len(unique_resolutions_minus_one) >= int(sys.argv[4]): # this is the case when the resolutions-1 are more than or equal to the number of genomes asked by the user (sys.argv[4])
                        for item in unique_resolutions_minus_one: #for each unique resolution-1
                            for key, value in acc_taxonomy.items(): #find it in acc_taxonomy
                                y = value.split(';')
                                if sys.argv[2] in y and item in y: #if both the <tax_level> and the resolution-1 are in the taxonomy string
                                    unique_resolution_minus_one_accessions.append(key) #add the accession to unique_resolution_minus_one_accessions, so it will contain all accessions of said resolution-1

                            for accession in unique_resolution_minus_one_accessions: #for each of these accessions
                                for key6,value6 in acc_score.items(): #find its score
                                    if accession == key6: #match the accession
                                        unique_resolution_minus_one_scores.append(value6) #add the score to a list for that resolution-1

                            for key7, value7 in acc_score.items():
                                if max(unique_resolution_minus_one_scores) == value7 and key7 in unique_resolution_minus_one_accessions:
                                    sub_acc_score.update( {key7 : value7} ) #make a sub-dictionary for failsafe, to match score with accession, only if accession is in res-1 list



                            unique_resolution_minus_one_scores.clear() #Clear the lists for the resolution-1
                            unique_resolution_minus_one_accessions.clear()
                        sorted_sub_acc_score = sorted(sub_acc_score.items(), key = lambda kv: kv[1], reverse=True) # sort the dictionary by value(score), but make it as tuple
                        sorted_dict = dict(sorted_sub_acc_score) #make the tuple dictionary again
                        for acc in list(sorted_dict.keys())[:int(sys.argv[4])]:#print <number> first keys from the sorted values dict
                            assemblies.write(acc[3:-2]+'\n')
                        sub_acc_score.clear()     # clear the sub-dict to proceed to next res-1



                    else: # this is the case when the resolution-1 are less than the number of genomes asked by the user (sys.argv[4])
                        f=0 # helps us determine number of loops
                        while(f<int(sys.argv[4])) and set(unique_resolutions_minus_one) != set(empty_res_minus_one): #The loop will run while fewer accessions than <number> have been picked and not all accession-1 have been emptied
                            for item in unique_resolutions_minus_one: #for each unique resolution-1
                                #print(item)
                                #if f < int(sys.argv[4]): # you NEED this extra fail-safe because while loop will run the whole for loop, so you have to check in every part of the for loop (so the for loop won't run once more before satisfying the condition of the while loop)
                                for key, value in acc_taxonomy.items(): #find it in acc_taxonomy
                                    y = value.split(';')
                                    if sys.argv[2] in y and item in y: #if both the <tax_level> and the resolution-1 are in the taxonomy string
                                        if key not in already_printed_accessions: #if the accession hasn't been printed already
                                            unique_resolution_minus_one_accessions.append(key) #add it to a list of unique accessions
                                for accession in unique_resolution_minus_one_accessions: #for each of these accessions
                                    for key6,value6 in acc_score.items():  #find its score
                                        if accession == key6: #match the accession
                                            unique_resolution_minus_one_scores.append(value6) #add the score to a list for that resolution-1



                                for key2, value2 in acc_score.items():
                                    if len(unique_resolution_minus_one_scores) != 0: # which means that not every accession of the resolution-1 has been printed
                                        if max(unique_resolution_minus_one_scores) == value2 and key2 in unique_resolution_minus_one_accessions: #find the element with that score, #failsafe for multiplicate scores, not to pick random accession with same score, but pick accession that exists in the unique resolution accessions
                                            sub_acc_score.update( {key2 : value2} ) #make a sub-dictionary for failsafe, to match score with accession, only if accession is in res-1 list
                                            break #break the loop so as not to append to the dict another accession with the same score (inside the same res-1)

                                    else: #if there are no more elements in that resolution-1
                                        if item not in empty_res_minus_one: #if it's not already in empty_res_minus_one
                                            empty_res_minus_one.append(item) #add it

                                unique_resolution_minus_one_accessions.clear() #clear the unique resolution-1 lists
                                unique_resolution_minus_one_scores.clear()

                                sorted_sub_acc_score = sorted(sub_acc_score.items(), key = lambda kv: kv[1], reverse=True) # sort the dictionary by value(score), but make it as tuple
                                sorted_dict = dict(sorted_sub_acc_score) #make the tuple dictionary again
                                for acc in list(sorted_dict.keys()):#print keys from the sorted values dict
                                    if f < int(sys.argv[4]):
                                        assemblies.write(acc[3:-2]+'\n')
                                        f=f+1
                                        already_printed_accessions.append(acc) #add the accession to already printed
                                sub_acc_score.clear()



                already_printed_accessions.clear()
                unique_resolutions_minus_one.clear()
                empty_res_minus_one.clear()
assemblies.close()



# This part is to create the assembliesnames text file for the contigs2orfs script
create_assembliesnames = str('while IFS= read -r line; do perl -n -e \'BEGIN { $matcher = pop } print if m/$matcher/g\' ' + sys.argv[1] + ' $line | perl -p -e \'s/^.*?_(.*?)\\.\\d+\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t.*?\\t(.*?)\\t.*/$1\\t$2/g\' >> ' + file_stem + '.assembliesnames ; done < ' + file_stem + '.assemblies')
if os.WEXITSTATUS(os.system(create_assembliesnames)) == 1:
    print('Error when writing the assembliesnames file. Exiting.')
    sys.exit(1)

#print('Genomes picked.')
assemblies_names.close()

# replace ';' in assembliesnames with '|'
with open(file_stem+".assembliesnames", "r+") as file:
    content = file.read()
    file.seek(0)
    file.write(content.replace(';', '|'))
#    file.truncate()

print('Converting RefSeq accesions to Genbank.')
rename1 = str("perl -p -i -e 's/^GCF_/GCA_/g' " + file_stem + ".assemblies")
#os.system(rename1)
if os.WEXITSTATUS(os.system(rename1)) == 1:
    print('Error when converting Refseq accessions to Genbank in assemblies file. Exiting.')
    sys.exit(1)

rename2 = str("perl -p -i -e 's/^GCF_/GCA_/g' " + file_stem + ".assembliesnames")
#os.system(rename2)
if os.WEXITSTATUS(os.system(rename2)) == 1:
    print('Error when converting Refseq accessions to Genbank in assembliesnames file. Exiting.')
    sys.exit(1)

print('All done!')
