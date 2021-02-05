#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.0
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Filter ancestral blocks (as if they were ancestral genes), keeping the largest blocks that account for minProp % of the genome
"""

import sys
import itertools

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("speciesTree",file), ("target",str), ("IN.ancGenesFiles",str), ("OUT.ancGenesFiles",str), ("minProp",int)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

for anc in targets:

		ls = []
		ll = []
		f = utils.myFile.openFile(arguments["IN.ancGenesFiles"] % phylTree.fileName[anc], "r")
		for l in f:
			ls.append(int(l.split("\t")[1]))
			ll.append(l)
		f.close()

		cutoff = utils.myMaths.myStats.getValueNX(sorted(ls), arguments["minProp"])
		f = utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc], "w")
		for (s,l) in itertools.izip(ls, ll):
			if s >= cutoff:
				print >> f, l,
			else:
				print >> f, "%s\t0\t\t\t" % anc
		f.close()

