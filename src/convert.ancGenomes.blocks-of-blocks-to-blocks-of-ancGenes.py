#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.0
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Convert a genome written as blocks of contigs to blocks of ancGenes
"""

import sys

import utils.myGraph
import utils.myFile
import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [("scaffoldsFile",file), ("contigsFile",file)], [], __doc__)

(diags,singletons) = utils.myGraph.loadIntegr(arguments["scaffoldsFile"])

ref = {}
f = utils.myFile.openFile(arguments["contigsFile"], "r")
for (i,l) in enumerate(f):
	ref[i+1] = l
f.close()

for (chrom,weights) in diags:
	li = []
	ls = []
	lw = []
	n = 0
	for (i,(c,s)) in enumerate(chrom):
		t = ref.pop(c)[:-1].split("\t")
		if i >= 1:
			lw.append(weights[i-1])
		n += len(t[2].split())
		if s > 0:
			li.append(t[2])
			ls.append(t[3])
			lw.append(t[4])
		else:
			li.extend(reversed(t[2].split()))
			ls.extend(-int(x) for x in reversed(t[3].split()))
			lw.extend(reversed(t[4].split()))
		
	print utils.myFile.myTSV.printLine([t[0], n, utils.myFile.myTSV.printLine(li, delim=" "), utils.myFile.myTSV.printLine(ls, delim=" "), utils.myFile.myTSV.printLine(lw, delim=" ")])

for c in singletons:
	print ref.pop(c),

# S'assure que tous les contigs ont ete employes
assert len(ref) == 0

