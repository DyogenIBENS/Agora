#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# AGORA v3.1
# python 3.5
# Copyright © 2006-2022 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021-2022 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Extract block adjacencies that are conserved in a given species set

    Usage:
        src/buildSynteny.pairwise-conservedAdjacencies.py example/data/Species.nwk A0 \
                -genesFiles=example/data/genes/genes.%s.list.bz2 \
                -ancGenesFiles=example/results/filtBlocks/best-pass1-all/blocks.%s.list.bz2 \
                -iniAncGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2 \
                -OUT.pairwise=example/results/pairwise/adjacencies-best-pass1-all/%s.list.bz2 \
                -LOG.pairwise=example/results/pairwise/adjacencies-best-pass1-all/%s.log.bz2
"""

import collections
import itertools
import multiprocessing
import sys
import time

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myGraph
from utils.myTools import file

# Arguments
arguments = utils.myTools.checkArgs(
	[("speciesTree",file), ("target",str)], \
	[("extantSpeciesFilter",str,""), \
	 ("genesFiles",str,""), ("ancGenesFiles",str,""), ("iniAncGenesFiles",str,""), ("OUT.pairwise",str,""),
	 ("anchorSize",int,2),
	 ("nbThreads", int, 0),
	 ("LOG.pairwise", str, "")],
	__doc__
)


def rewriteGenome(genome, dic):
	newGenome = {}
	for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
		tmp = [(dic[(chrom,i)],gene.strand) for (i,gene) in enumerate(genome.lstGenes[chrom]) if (chrom,i) in dic]
		tmp = [(i,s*sa) for ((i,sa),s) in tmp]
		if len(tmp) > 0:
			newGenome[chrom] = [i for (i,l) in itertools.groupby(tmp)]
	return newGenome

def getExtremities(genome):
	extr1 = {}
	extr2 = {}
	for (chrom,l) in genome.items():
		(i0,s0) = l[0]
		(i1,s1) = l[-1]
		extr1[(i0,s0)] = (chrom,1)
		extr2[(i1,s1)] = (chrom,1)
		extr1[(i1,-s1)] = (chrom,-1)
		extr2[(i0,-s0)] = (chrom,-1)
	return (extr1, extr2)

def getAllAdj(anc):
	allAdj = collections.defaultdict(list)
	anchorSize = arguments["anchorSize"]
	for x in dicGenomes[anc].lstGenes.values():
		if (len(x) >= 2) and (len(x) < anchorSize):
			anchorSize = len(x)

	log = arguments["LOG.pairwise"] % phylTree.fileName[anc]
	f = utils.myFile.openFile(log, "w")
	for esp in sorted(listSpecies):

		dicA = {}
		dicM = {}
		stats = []

		for (n,((c1,d1),(c2,d2),da)) in enumerate(utils.myGraph.calcDiags(dicGenomes[esp], dicGenomes[anc], genesAnc[phylTree.dicParents[anc][esp]], orthosFilter=utils.myGraph.OrthosFilterType.InBothSpecies, minChromLength=anchorSize)):
			if len(da) < anchorSize:
				continue
			print("DIAG", anc, esp, n, (c1,c2), len(da), (d1,d2,da), file=f)
			for ((i1,s1),(i2,s2)) in zip(d1, d2):
				dicM[(c1,i1)] = (n,s1)
				dicA[(c2,i2)] = (n,s1)
			stats.append(len(da))

		newGA = rewriteGenome(dicGenomes[anc], dicA)
		# The blocs selected so far
		notdup = set()
		for cA in newGA:
			notdup.update(x[0] for x in newGA[cA])
			print(anc, esp, "ANC", cA, len(newGA[cA]), newGA[cA], file=f)
		newGM = rewriteGenome(dicGenomes[esp], dicM)

		(extr1,extr2) = getExtremities(newGA)

		na = 0
		for (cM,l) in newGM.items():
			print(anc, esp, "MOD", cM, len(newGM[cM]), newGM[cM], file=f)
			# In case there is a segmental duplication, choose the same block as in the ancestor
			l = [x for x in l if x[0] in notdup]
			print(anc, esp, "FMOD", cM, len(l), l, file=f)
			for (x1,x2) in utils.myTools.myIterator.slidingTuple(l):
				if (x1 in extr2) and (x2 in extr1):
					(i1,s1) = extr2[x1]
					(i2,s2) = extr1[x2]
					if i1 == i2:
						print("LOOP", extr2[x1], extr1[x2], file=f)
						continue
					print("ADJ", anc, esp, (x1,x2), (extr2[x1],extr1[x2]), file=f)
					na += 1
					if i1 < i2:
						allAdj[ ((i1,s1),(i2,s2)) ].append(esp)
					else:
						allAdj[ ((i2,-s2),(i1,-s1)) ].append(esp)
		print("Gene order comparison between %s and %s ..." % (anc,esp), utils.myMaths.myStats.txtSummary(stats), "%d adjacencies / %d blocks" % (na, len(newGA)), "(anchor size: %d)" % anchorSize, file=sys.stderr)
	f.close()

	# -1 is the outgroup species, 1,2,3... are the descendants
	ind = dict.fromkeys(set(phylTree.outgroupSpecies[anc]).intersection(listSpecies), -1)
	for (i,(x,_)) in enumerate(phylTree.items[anc]):
		ind.update(dict.fromkeys(set(phylTree.species[x]).intersection(listSpecies), i+1))

	res = arguments["OUT.pairwise"] % phylTree.fileName[anc]
	f = utils.myFile.openFile(res, "w")
	for ancPair in allAdj:

		# Compute the weight (number of comparisons that support this adjacency)
		weights = collections.defaultdict(int)
		for esp in allAdj[ancPair]:
			weights[ind[esp]] += 1
		weight = sum(x*y for (x,y) in itertools.combinations(list(weights.values()), 2))

		if weight == 0:
			# None of the adjacencies is conserved through this ancestor
			continue

		# to debug the weight: "|".join("%s/%d" % (esp,ind[esp]) for esp in sorted(allAdj[ancPair]))
		print(utils.myFile.myTSV.printLine(
			list(ancPair[0] + ancPair[1]) +
			[weight]
		), file=f)
	f.close()



phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

(listSpecies, targets, accessoryAncestors) = phylTree.getTargetsForPairwise(arguments["target"], arguments["extantSpeciesFilter"])

dicGenomes = {}
for e in sorted(listSpecies):
	dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[e])

genesAnc = {}
for anc in sorted(targets.union(accessoryAncestors)):
	genesAnc[anc] = utils.myGenomes.Genome(arguments["iniAncGenesFiles"] % phylTree.fileName[anc])
for anc in sorted(targets):
	dicGenomes[anc] = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc], ancGenes=genesAnc[anc], withDict=False)

toStudy = collections.defaultdict(list)
for (e1,e2) in itertools.combinations(listSpecies, 2):
	for anc in targets.intersection(phylTree.dicLinks[e1][e2][1:-1]):
		toStudy[anc].append((e1,e2))

start = time.time()
n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
multiprocessing.Pool(n_cpu).map(getAllAdj, sorted(targets))
print("Elapsed time:", (time.time() - start), file=sys.stderr)

