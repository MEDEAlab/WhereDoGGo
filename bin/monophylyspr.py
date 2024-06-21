#!/usr/bin/env python

#The MIT License (MIT) Copyright (c) 2024 George Kolyfetis & Panagiotis Adam
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#Function
#This script checks if leaves from a given list are monophyletic in a phylogoeny. If they are, the branch is pruned and regradted to all possible positions in the phylogeny.

#NOTE 1: All code was written and tested on Intel macOS and Ubuntu. Please report any issues.

#Dependencies
#1) ETE3 (https://anaconda.org/etetoolkit/ete3 or https://pypi.org/project/ete3/)

import os
import sys

#Check if required non-standard libraries are installed.
import importlib.util
nonstandardlibraries = {"ete3": "https://anaconda.org/etetoolkit/ete3 or https://pypi.org/project/ete3/"}
for nstlobject,link in nonstandardlibraries.items():
    if importlib.util.find_spec(nstlobject) is not None:
        pass
    else:
        print('Library ' + nstlobject + ' not installed. Download it from: ' + link + '. Exiting.')
        sys.exit(1)

from ete3 import Tree

print('#Script: monophylyspr.py')
print('#Version: v20240621')
print('#Usage: python monophylyspr.py <input_tree> <cluster>')
print('#<input_tree> must be the input phylogenetic tree in Newick format. Branch supports need to be single values. (required)')
print('#<cluster> must be a text file containing leaf names (one per line). These leaves are checked if they are monophyletic in <input_tree>. If they are, their branch is pruned and regrafted in all possible positions on the tree. (required)')
print('#For more information refer to the comments in the script and/or the Github page.')

# Check if the correct number of arguments is given
if len(sys.argv) == 3:
    print ('Two arguments found. Proceeding.')
else:
    print('Wrong number of arguments given. Exiting.')
    sys.exit(1)

# checkpoint for input tree file
check_file = os.path.isfile(sys.argv[1])
if check_file == True:
    print('Input tree file found. Proceeding.')
else:
    print('Input tree file not found. Exiting.')
    sys.exit(1)
#TODO: Check for Newick format and single branch supports?

# checkpoint for leaf names file
check_file = os.path.isfile(sys.argv[2])
if check_file == True:
    print('File with leaf names found. Proceeding.')
else:
    print('File with leaf names not found. Exiting.')
    sys.exit(1)
#TODO: Check for formatting?

def remove_identical_trees(input_tree_path, multi_newick_path, output_path):
    input_tree = Tree(input_tree_path, format=0)
    with open(multi_newick_path, 'r') as file:
        multi_newick_trees = [line.strip() for line in file]
    alreadyfound = False
    for tree_str in multi_newick_trees:
        checktree = Tree(tree_str, format=0)
        rf_distance = str(input_tree.robinson_foulds(checktree, unrooted_trees=True)[0])
        if rf_distance == "0" and alreadyfound == False:
            alreadyfound = True
        elif rf_distance == "0" and alreadyfound == True:
            multi_newick_trees.remove(tree_str)
    with open(output_path, 'w') as file:
        file.write('\n'.join(multi_newick_trees))

#Remove files from previous runs
print('Removing files with names identical to the output.')
removal = ('rm -r temp_tree.tree SPR*.tree* 2> /dev/null')
os.system(removal)

input_cluster=list()

for i in open(sys.argv[2], 'r'): #this is a text file with the leaf names to be checked, pruned, and regrafted, one per line
    input_cluster.append(i.strip())

g = open('temp_tree.tree', 'a') #we need to cp the tree to another file so as to use it a second time for the second pruned tree

with open(sys.argv[1], 'r') as tree:
    for line in tree:
        g.write(line)

g.close()

t = Tree(sys.argv[1]) #this tree needs to have single branch support values
o = Tree('temp_tree.tree', format=0)

#create a list with all tree leaves except the input taxa
leaves_except_input = list()
for leaf in t:
    if leaf.name not in input_cluster:
        leaves_except_input.append(leaf.name)

print('Checking monophyly of the list of leaves.')
# check the monophyly of the cluster given
if t.check_monophyly(values=input_cluster, target_attr="name") == (True, 'monophyletic', set()):
    o.prune(input_cluster, preserve_branch_length=True) #this is the input cluster
    print('Given leaves are monophyletic. Proceeding with all possible pruning and regrafting moves.')
    y=0
    t.prune(leaves_except_input) #tree(t) is now missing the input cluster and is ready to be modified multiple times
    l=0
    d=0
## EXTERNAL
    for node in t.traverse():
        l=l+1
        if node.name == '':
            if node.is_root():
                node.add_feature('number', -2)
            else:
                node.add_feature('number', y)
                y=y+1
        else:
            node.add_feature('number', -1)
    #print('Number of nodes (internal + leaves):' + ' ' + str(l) + '\n')

    for node in t.traverse('postorder'):
        t = Tree(sys.argv[1])
        t.prune(leaves_except_input)
        o = Tree('temp_tree.tree')
        o.prune(input_cluster, preserve_branch_length=True)

        if node.number == -1: #aka is external
            d = d+1
            A = t.search_nodes(name=node.name)[0]
            A.add_child(o)
            A.add_child(name=node.name)
            t.write(outfile='SPR_tree%s.tree' % d)
## INTERNAL AND ROOT
    for node in t.traverse():
        l=l+1
        if node.name == '':
            if node.is_root():
                node.add_feature('number', -2)
            else:
                node.add_feature('number', y)
                y=y+1
        else:
            node.add_feature('number', -1)

    for node in t.traverse('postorder'):
        t.prune(leaves_except_input)
        o = Tree('temp_tree.tree')
        o.prune(input_cluster, preserve_branch_length=True)

        if node.number != -1 and node.number != -2: #aka is internal and not root
            d = d+1
            G = node.populate(1, names_library='G')
            for node in t.traverse():
                if node.name == 'G':
                    node.add_child(o)
                    node.delete()
                    #t.convert_to_ultrametric(tree_length=1)
            t.write(outfile='SPR_tree%s.tree' % d)

        elif node.number == -2: #aka is root
            d = d+1
            b = Tree()
            b.add_child(node, dist = 0.05)
            b.add_child(o)
            b.write(outfile='SPR_tree%s.tree' % d)

    #print('Number of created SPR trees:' + ' ' + str(d) + '\n')

elif len(input_cluster) == 1:
    o.prune(input_cluster, preserve_branch_length=True) #this is the input cluster
    print('Given leaves are monophyletic. Proceeding with all possible pruning and regrafting moves.')
    y=0
    t.prune(leaves_except_input) #tree(t) is now missing the input cluster and is ready to be modified multiple times
    l=0

    d=0
## EXTERNAL
    for node in t.traverse():
        l=l+1
        if node.name == '':
            if node.is_root():
                node.add_feature('number', -2)
            else:
                node.add_feature('number', y)
                y=y+1
        else:
            node.add_feature('number', -1)
    #print('Number of nodes (internal + leaves):' + ' ' + str(l) + '\n')

    for node in t.traverse('postorder'):
        c = Tree()
        t.prune(leaves_except_input)
        o = Tree('temp_tree.tree')
        o.prune(input_cluster, preserve_branch_length=True)

        if node.number == -1: #aka is external
            d = d+1
            c.add_child(name=node.name, dist = 0.05)
            c.add_child(name=input_cluster[0])
            node.add_child(c)
            node.delete()
            t.write(outfile='SPR_tree%s.tree' % d)

## INTERNAL AND ROOT
    for node in t.traverse():
        l=l+1
        if node.name == '':
            if node.is_root():
                node.add_feature('number', -2)
            else:
                node.add_feature('number', y)
                y=y+1
        else:
            node.add_feature('number', -1)
    for node in t.traverse('postorder'):
        t.prune(leaves_except_input)
        o = Tree('temp_tree.tree')
        o.prune(input_cluster, preserve_branch_length=True)

        if node.number != -1 and node.number != -2: #aka is internal and not root
            d = d+1
            G = node.populate(1, names_library='G')
            for node in t.traverse():
                if node.name == 'G':
                    node.add_child(name=input_cluster[0])
                    node.delete()
                    #t.convert_to_ultrametric(tree_length=1)
            t.write(outfile='SPR_tree%s.tree' % d)
        elif node.number == -2: #aka is root
            d = d+1
            e = Tree()
            e.add_child(node, dist = 0.05)
            e.add_child(name=input_cluster[0])
            e.write(outfile='SPR_tree%s.tree' % d)

    #print('Number of created SPR trees:' + ' ' + str(d) + '\n')

else:
    print('Given leaves are not monophyletic. Cannot perform pruning and regrafting.')

print('Removing temporary files and placing all trees in a single file. Any branch supports will also be removed.')
os.remove("temp_tree.tree")
os.system('cat ' + sys.argv[1] + ' >> SPR_trees.treels && for file in SPR_tree*.tree ; do cat "$file" >> SPR_trees.treels && echo "" >> SPR_trees.treels && rm "$file" ; done && perl -p -i -e \'s/^\\n//g\' SPR_trees.treels') # the one-liner removes empty lines. #&& perl -p -i -e \'s/\\)\\d+\\:/\\)\\:/g\' SPR_trees.treels

print('Removing duplicate trees.')
remove_identical_trees(sys.argv[1], 'SPR_trees.treels', 'SPR_trees.treels')

print('All done!')
