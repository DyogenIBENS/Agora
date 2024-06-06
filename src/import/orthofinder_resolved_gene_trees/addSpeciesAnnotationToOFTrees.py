#!/bin/env python

from ete3 import Tree
import sys
import os
import argparse
__doc__ = """
script name: addSpeciesAnnotationToOFTrees.py

Example : python ./addSpeciesAnnotationToOFTrees.py -of_trees OrthoFinder/Results_Apr07/Resolved_Gene_Trees  -outdir Agora_data/Resolved_Gene_trees



"""
if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    PARSER.add_argument('-of_trees', '--orthofinder_trees', help='Path to orthofinder2 Resolved gene trees directory',\
                         required=True)
    PARSER.add_argument('-outdir', '--outdir', help='Path to AGORA orthoGroups directory',\
                         required=True)

    ARGS = vars(PARSER.parse_args())

dirTrees= ARGS['orthofinder_trees']
outOG=ARGS['outdir']
print("Parsing:", dirTrees, file=sys.stderr)

for filename in os.listdir(dirTrees):
    f = os.path.join(dirTrees, filename)
    base = os.path.splitext(filename)[0]
    filename2=base+'.nhx'
    of = os.path.join(outOG,filename2)
    
    tree_file=Tree(f,format=1)

    for n in tree_file.traverse():
        if n.is_leaf():
            removal = "_"        
            reverse_removal = removal[::-1]
            replacement = " "
            reverse_replacement = replacement[::-1]
            newstr = n.name[::-1].replace(reverse_removal, reverse_replacement, 1)[::-1]
            tmp=newstr.split()
            n.name=tmp[1]
            n.S=tmp[0]

    outfile=open(of,"w") 
    print(tree_file.write(format=9,  features=["S"]), file=outfile)
    outfile.close()
    print(f+" done... see "+of,file=sys.stderr)
