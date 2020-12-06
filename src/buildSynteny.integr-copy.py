#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.0
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Just a script to copy the ancestral blocks of robusts genes (skeleton) in the good directory to continue the procedure
"""

import time
import sys

import utils.myFile
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs([("speciesTree", file), ("target", str)],
                                    [("IN.ancBlocks", str, ""), ("OUT.ancBlocks", str, "")],
                                    __doc__
                                    )

start = time.time()
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

for anc in targets:
    fi = utils.myFile.openFile(arguments["IN.ancBlocks"] % phylTree.fileName[anc], "r")
    fo = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    for l in fi:
        print >> fo, l,
    fo.close()
    fi.close()

print >> sys.stderr, "Time elapsed:", time.time() - start
