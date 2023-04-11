#!/usr/bin/env python3

import sys
import argparse
import os
import re
__doc__ = """
script name: convert_hogs_sp.py

Example : python ./convert_hogs_sp.py -of_hogs OrthoFinder/Results_Apr07/Phylogenetic_Hierarchical_Orthogroups  -outdir Agora_data/orthoGroups



"""


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    PARSER.add_argument('-of_hogs', '--orthofinder_hogs', help='Path to orthofinder2 Hogs directory',\
                         required=True)
    PARSER.add_argument('-outdir', '--outdir', help='Path to AGORA orthoGroups directory',\
                         required=True)

    ARGS = vars(PARSER.parse_args())


dirHogs= ARGS['orthofinder_hogs']
outOG=ARGS['outdir']
print("Parsing:", dirHogs, file=sys.stderr)


for filename in os.listdir(dirHogs):
    f = os.path.join(dirHogs, filename)
    base = os.path.splitext(filename)[0]
    filename2="orthologyGroups."+base+'.list'
    of = os.path.join(outOG,filename2)
    # checking if it is a file
    if os.path.isfile(f) :
        Lines=open(f,"r").read().splitlines()
        if Lines[0].split("\t")[0] == "HOG":
            print("formating:",f," for AGORA",file=sys.stderr)
            outfile=open(of,"w") 
            for l in Lines[1:]:
                    og=re.split("[\t,]+",l)[3:]
                    print(*og, sep = " ", file=outfile)
            outfile.close()
            print("orthoGroup file for AGORA: ", of, file=sys.stderr) 
