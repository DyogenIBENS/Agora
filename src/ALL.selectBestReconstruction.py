#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	For each ancestor, select the reconstruction that has the highest N50
"""

import sys
import itertools

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs(
        [("speciesTree",file), ("target",str), ("IN.ancBlocks", utils.myTools.ParamList(str, 1))],
        [("OUT.ancBlocks", str, "")],
        __doc__
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

for anc in targets:

    stats = {}
    lines = {}
    for ancBlocksFileName in arguments["IN.ancBlocks"]:
        ls = []
        ll = []
        f = utils.myFile.openFile(ancBlocksFileName % phylTree.fileName[anc], "r")
        for l in f:
            ls.append(int(l.split("\t")[1]))
            ll.append(l)
        f.close()
        ls.sort()
        lines[ancBlocksFileName] = ll
        stats[ancBlocksFileName] = utils.myMaths.myStats.getValueNX(ls, 50)
        print >> sys.stderr, "%s @ %s: %s" % (anc, ancBlocksFileName, utils.myMaths.myStats.txtSummary(ls))

    bestBlocks = max(arguments["IN.ancBlocks"], key=stats.get)
    print >> sys.stderr, "-> Best for %s is %s" % (anc, bestBlocks)

    fo = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    for l in lines[bestBlocks]:
        print >> fo, l,
    fo.close()

