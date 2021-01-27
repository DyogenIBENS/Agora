#! /usr/bin/env python

__doc__ = """
	Filtre les contigs d'une taille minimale
"""

import sys

import utils.myFile
import utils.myTools
import utils.myPhylTree

arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("target",str), ("IN.ancGenesFiles",str), ("OUT.ancGenesFiles",str), ("minLength",int)], [], __doc__ )

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
target = phylTree.officialName[arguments["target"]]

for anc in phylTree.listAncestr:
	if phylTree.dicParents[anc][target] == target:
		fi = utils.myFile.openFile(arguments["IN.ancGenesFiles"] % phylTree.fileName[anc], "r")
		fo = utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc], "w")
		for l in fi:
			if int(l.split("\t")[1]) >= arguments["minLength"]:
				print >> fo, l,
			else:
				print >> fo, "%s\t0\t\t\t" % anc
		fi.close()
		fo.close()

