#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.0
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Simple script to copy ancestral blocks to another directory. Useful to combine
    different reconstructions (on different ancestors) in one directory

    Usage:
        src/buildSynteny.integr-copy.py example/data/Species.nwk A0 \
                -IN.ancBlocks=example/results/ancBlocks/best-pass1/blocks.%s.list.bz2 \
                -OUT.ancBlocks=example/results/filtBlocks/best-pass1-all/blocks.%s.list.bz2
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
