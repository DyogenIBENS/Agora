#! /usr/bin/env python
#  -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """ Convert a species phylogenetic tree from/to newick or phyltree text format

     usage:
        ALL.convertNewickTree.py Species.conf -fromNewick > Species.nwk
        ALL.convertNewickTree.py Species.nwk +fromNewick > Species.conf """

import utils.myFile as myFile
import utils.myPhylTree as myPhylTree
import utils.myTools as myTools


arguments = myTools.checkArgs([("phylTree.conf",file)], [("fromNewick",bool,True)], __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])


if arguments["fromNewick"]:

    # print into phylTree format, with tabulations
    def do(node, indent):
        node = node.replace("*", "")
        node = node.replace(".", " ")
        names = myFile.myTSV.printLine([node] + [x for x in phylTree.commonNames.get(node,"") if isinstance(x, str) and (x != node)], delim="|")
        if node.replace(" ",".") in phylTree.listSpecies :
            print ("\t" * indent) + str(names)
        elif node in phylTree.items:
            print ("\t" * indent) + str(names) + "\t" + str(int(phylTree.ages[node]))
            for (f,_) in phylTree.items[node]:
                do(f, indent+1)
    do(phylTree.root, 0)

else:
    # return the tree into the newick tree
    def convertToFlatFile(anc):
        a = phylTree.fileName[anc] # anc.replace(' ', '.')
        if anc in phylTree.listSpecies:
            return a
        else:
            return "(" + ",".join([convertToFlatFile(e) + ":" + str(l) for (e,l) in phylTree.items[anc]]) + ")%s" % a #" |%d" % (a,phylTree.ages[anc])
    print convertToFlatFile(phylTree.root), ";"
