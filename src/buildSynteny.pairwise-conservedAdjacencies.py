#! /usr/bin/env python

__doc__ = """
	Extrait les adjacences de contigs conservees
"""

import sys
import itertools
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myDiags

# Arguments
arguments = utils.myTools.checkArgs(
	[("phylTree.conf",file), ("target",str)], \
	[("genesFiles",str,""), ("ancGenesFiles",str,""),
	("iniAncGenesFiles",str,""), ("usedSpecies",str,""), ("anchorSize",int,2), ("verbose",bool,True)],
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
	for (chrom,l) in genome.iteritems():
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
	for x in dicGenomes[anc].lstGenes.itervalues():
		if (len(x) >= 2) and (len(x) < anchorSize):
			anchorSize = len(x)

	for esp in listSpecies:

		dicA = {}
		dicM = {}
		stats = []

		for (n,((c1,d1),(c2,d2),da)) in enumerate(utils.myDiags.calcDiags(dicGenomes[esp], dicGenomes[anc], genesAnc[phylTree.dicParents[anc][esp]], orthosFilter=utils.myDiags.OrthosFilterType.InBothSpecies, minChromLength=anchorSize)):
			if len(da) < anchorSize:
				continue
			if arguments["verbose"]:
				print >> sys.stderr, "DIAG", anc, esp, n, (c1,c2), len(da), (d1,d2,da)
			for ((i1,s1),(i2,s2)) in itertools.izip(d1, d2):
				dicM[(c1,i1)] = (n,s1)
				dicA[(c2,i2)] = (n,s1)
			stats.append(len(da))

		newGA = rewriteGenome(dicGenomes[anc], dicA)
		# Liste des blocs choisis chez l'ancetre
		notdup = set()
		for cA in newGA:
			notdup.update(x[0] for x in newGA[cA])
			if arguments["verbose"]:
				print >> sys.stderr, anc, esp, "ANC", cA, len(newGA[cA]), newGA[cA]
		newGM = rewriteGenome(dicGenomes[esp], dicM)

		(extr1,extr2) = getExtremities(newGA)

		na = 0
		for (cM,l) in newGM.iteritems():
			if arguments["verbose"]:
				print >> sys.stderr, anc, esp, "MOD", cM, len(newGM[cM]), newGM[cM]
			# Permet de selectionner lors d'une duplication segmentale, le meme bloc que chez l'ancetre
			l = [x for x in l if x[0] in notdup]
			if arguments["verbose"]:
				print >> sys.stderr, anc, esp, "FMOD", cM, len(l), l
			for (x1,x2) in utils.myTools.myIterator.slidingTuple(l):
				if (x1 in extr2) and (x2 in extr1):
					(i1,s1) = extr2[x1]
					(i2,s2) = extr1[x2]
					if i1 == i2:
						if arguments["verbose"]:
							print >> sys.stderr, "LOOP", extr2[x1], extr1[x2]
						continue
					if arguments["verbose"]:
						print >> sys.stderr, "ADJ", anc, esp, (x1,x2), (extr2[x1],extr1[x2])
					na += 1
					if i1 < i2:
						allAdj[ ((i1,s1),(i2,s2)) ].append(esp)
					else:
						allAdj[ ((i2,-s2),(i1,-s1)) ].append(esp)
		print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (anc,esp), utils.myMaths.myStats.txtSummary(stats), "%d adjacences / %d blocs" % (na, len(newGA)), "(ancre: %d)" % anchorSize

	# -1 designe l'outgroup, 1,2,3... designent les descendants
	ind = dict.fromkeys(set(phylTree.outgroupSpecies[anc]).intersection(listSpecies), -1)
	for (i,(x,_)) in enumerate(phylTree.items[anc]):
		ind.update(dict.fromkeys(set(phylTree.species[x]).intersection(listSpecies), i+1))

	for ancPair in allAdj:

		# Calcul des poids
		weights = collections.defaultdict(int)
		for esp in allAdj[ancPair]:
			weights[ind[esp]] += 1
		weight = sum(x*y for (x,y) in itertools.combinations(weights.values(), 2))

		print utils.myFile.myTSV.printLine(
			[anc] + list(ancPair[0] + ancPair[1]) +
			[None, None, weight,
			"|".join("%s/%d" % (esp,ind[esp]) for esp in sorted(allAdj[ancPair]))]
		)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

targets = phylTree.getTargetsAnc(arguments["target"])
listSpecies = phylTree.getTargetsSpec(arguments["usedSpecies"] if len(arguments["usedSpecies"]) > 0 else phylTree.root)

dicGenomes = {}
for e in listSpecies:
	dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[e])

genesAnc = {}
for anc in targets:
	genesAnc[anc] = utils.myGenomes.Genome(arguments["iniAncGenesFiles"] % phylTree.fileName[anc])
	dicGenomes[anc] = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc], ancGenes=genesAnc[anc], withDict=False)

for anc in [phylTree.dicParents[e][a] for (e,a) in itertools.product(listSpecies, targets)]:
	if anc not in genesAnc:
		genesAnc[anc] = utils.myGenomes.Genome(arguments["iniAncGenesFiles"] % phylTree.fileName[anc])

toStudy = collections.defaultdict(list)
for (e1,e2) in itertools.combinations(listSpecies, 2):
	for anc in targets.intersection(phylTree.dicLinks[e1][e2][1:-1]):
		toStudy[anc].append((e1,e2))

for anc in targets:
	getAllAdj(anc)


