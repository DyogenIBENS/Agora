#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Just a script to copy ancestral blocks to another directory to continue the procedure
"""

import time
import sys

import utils.myFile
import utils.myPhylTree
import utils.myTools
from utils.myTools import file

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
        fo.write(l)
    fo.close()
    fi.close()

print("Time elapsed:", time.time() - start, file=sys.stderr)
