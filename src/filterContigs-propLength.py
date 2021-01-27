#! /usr/bin/env python

__doc__ = """
	Filtre les contigs d'une taille minimale
"""

import sys
import itertools

import utils.myFile
import utils.myMaths
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("target",str), ("IN.ancGenesFiles",str), ("OUT.ancGenesFiles",str), ("minProp",int)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
target = phylTree.officialName[arguments["target"]]

for anc in phylTree.listAncestr:
	if phylTree.dicParents[anc][target] == target:

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

