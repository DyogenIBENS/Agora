#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.1
# python 3.5
# Copyright © 2006-2022 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021-2022 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Filter ancestral blocks (as if they were ancestral genes), keeping the largest blocks that account for minProp % of the genome

    Usage:
        src/ALL.filterBlocks-propLength.py example/data/Species.nwk A0 \
                example/results/filtBlocks/best-pass1-all/blocks.%s.list.bz2 \
                example/results/filtBlocks/best-pass1-propLength-50/blocks.%s.list.bz2 \
                50
"""

import sys
import itertools

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree
from utils.myTools import file

arguments = utils.myTools.checkArgs( [("speciesTree",file), ("target",str), ("IN.ancGenesFiles",str), ("OUT.ancGenesFiles",str), ("minProp",int)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

for anc in sorted(targets):

		ls = []
		ll = []
		f = utils.myFile.openFile(arguments["IN.ancGenesFiles"] % phylTree.fileName[anc], "r")
		for l in f:
			ls.append(int(l.split("\t")[1]))
			ll.append(l)
		f.close()

		cutoff = utils.myMaths.myStats.getValueNX(sorted(ls), arguments["minProp"])
		f = utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc], "w")
		for (s,l) in zip(ls, ll):
			if s >= cutoff:
				f.write(l)
			else:
				print("%s\t0\t\t\t" % anc, file=f)
		f.close()

