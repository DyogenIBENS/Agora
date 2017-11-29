#! /usr/bin/env python

__doc__ = """
	Lance walktrap et renvoie le meilleur clustering
"""

import sys
import operator
import collections

import utils.myTools
import utils.walktrap

arguments = utils.myTools.checkArgs( [("data",file)], [("randomWalksLength",int,4), ("qualityFunction",int,[2,1,3])], __doc__ )

print >> sys.stderr, "Loading ...",
edges = collections.defaultdict(dict)
f = utils.myFile.openFile(arguments["data"], "r")
for l in f:
	t = l.split()
	edges[t[0]][t[1]] = edges[t[1]][t[0]] = float(t[2])
print >> sys.stderr, "OK"
f.close()

res = utils.walktrap.doWalktrap(edges, showProgress=True, randomWalksLength=arguments["randomWalksLength"], qualityFunction=arguments["qualityFunction"])

print >> sys.stderr, "Printing ...",
for (nodes,cuts,_,dend) in res:
	if len(cuts) == 0:
		print " ".join(nodes)
	else:
		(alpha,score) = max(cuts, key=operator.itemgetter(1))
		(clust,lonely) = dend.cut(alpha)
		for l in clust:
			print " ".join(l)
		if len(lonely) > 0:
			print " ".join(lonely)
print >> sys.stderr, "OK"

