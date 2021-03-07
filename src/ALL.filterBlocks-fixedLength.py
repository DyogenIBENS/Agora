#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Filter ancestral blocks (as if they were ancestral genes), using a size threshold (absolute value)s
"""

import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree
from utils.myTools import file

arguments = utils.myTools.checkArgs( [("speciesTree",file), ("target",str), ("IN.ancGenesFiles",str), ("OUT.ancGenesFiles",str), ("minLength",int)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

for anc in sorted(targets):
		fi = utils.myFile.openFile(arguments["IN.ancGenesFiles"] % phylTree.fileName[anc], "r")
		fo = utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc], "w")
		for l in fi:
			if int(l.split("\t")[1]) >= arguments["minLength"]:
				fo.write(l)
			else:
				print("%s\t0\t\t\t" % anc, file=fo)
		fi.close()
		fo.close()

