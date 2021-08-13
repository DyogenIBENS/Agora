#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.0
# python 3.5
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    For each ancestor, select the reconstruction that has the highest G50

    Usage:
        src/ALL.selectBestReconstruction.py example/data/Species.nwk A0 \
                -OUT.ancBlocks=example/results/ancBlocks/best-pass1.best-pass2/blocks.%s.list.bz2 \
                example/results/ancBlocks/best-pass1.denovo-all.asAncGenes/blocks.%s.list.bz2 \
                example/results/ancBlocks/best-pass1.denovo-propLength-50.fillin-all.fusion-all.insertion-all.asAncGenes/blocks.%s.list.bz2 \
                example/results/ancBlocks/best-pass1.denovo-propLength-70.fillin-all.fusion-all.insertion-all.asAncGenes/blocks.%s.list.bz2 \
                example/results/ancBlocks/best-pass1.denovo-fixedLength-20.fillin-all.fusion-all.insertion-all.asAncGenes/blocks.%s.list.bz2 \
                example/results/ancBlocks/best-pass1.denovo-fixedLength-50.fillin-all.fusion-all.insertion-all.asAncGenes/blocks.%s.list.bz2
"""

import sys
import itertools

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree
from utils.myTools import file

arguments = utils.myTools.checkArgs(
        [("speciesTree",file), ("target",str), ("IN.ancBlocks", utils.myTools.ParamList(str, 1))],
        [("OUT.ancBlocks", str, "")],
        __doc__
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

for anc in sorted(targets):

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
        print("%s @ %s: %s" % (anc, ancBlocksFileName, utils.myMaths.myStats.txtSummary(ls)), file=sys.stderr)

    bestBlocks = max(arguments["IN.ancBlocks"], key=stats.get)
    print("-> Best for %s is %s" % (anc, bestBlocks), file=sys.stderr)

    fo = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    for l in lines[bestBlocks]:
        fo.write(l)
    fo.close()

