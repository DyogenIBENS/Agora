#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Extract all pairs of genes (adjacencies) conserved between each pair of species.
"""

import collections
import itertools
import sys
import time

import utils.myMaths
import utils.myPhylTree
import utils.myGenomes


# Arguments
arguments = utils.myTools.checkArgs(
	[("speciesTree",file), ("target",str)], \
	[("extantSpeciesFilter",str,""), ("genesFiles",str,""), ("ancGenesFiles",str,""), ("OUT.pairwise",str,"")],
	__doc__
)


st = start = time.time()
#Species tree
#############
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

# Species to use
################

(listSpecies, listAncestors, accessoryAncestors) = phylTree.getTargetsForPairwise(arguments["target"], arguments["extantSpeciesFilter"])

def revPair((g1, g2)):
	return ((g2[0],-g2[1]),(g1[0],-g1[1]))

dicAncMod = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(list)))
dicModAnc = collections.defaultdict(list)

# We override intern() in order to be able to clear its cache once all the loading is done
name_hash = {}
def myintern(s):
	return name_hash.setdefault(s, s)
utils.myGenomes.intern = myintern

genesAnc = {}
for anc in listAncestors.union(accessoryAncestors):
	ancGenes = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])
	genesAnc[anc] = {k: v.index for (k,v) in ancGenes.dicGenes.iteritems()}
	del ancGenes

print >> sys.stderr, "time for loading", time.time() - start
start = time.time()

todo = {}
for esp in listSpecies:
	lanc = []
	anc = esp
	while anc in phylTree.parent:
		(par,_) = phylTree.parent[anc]
		if par in genesAnc:
			lanc.append((par, genesAnc[par], dicAncMod[par][anc]))
		anc = par
	todo[esp] = lanc
del genesAnc

def extractPairsFromSpecies(esp):
	genome = utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[esp], withDict=False)

	print >> sys.stderr, "Extraction of pairs of genes from %s " % esp, "...",
	
	for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[utils.myGenomes.ContigType.Scaffold]:
		chrom = genome.lstGenes[chrom]
		
		if len(chrom) < 2:
			continue
		chrom = [(None, (gene.names[-1], gene.strand)) for gene in chrom]
		
		for (anc,dica,subdicAncMod) in todo[esp]:
			# Updating the chromosome under the new ancestor, the list keeps on shrinking
			chrom = [((dica.pop(x[0]), x[1]), x) for (_,x) in chrom if x[0] in dica]
			if len(chrom) < 2:
				break
			
			(ga1, gm1) = chrom[0]
			for (ga2, gm2) in itertools.islice(chrom, 1, None):
				if ga1[0] == ga2[0]:
					ga1 = ga2
					gm1 = gm2
					continue
				# We only keep the pair in the right direction
				if ga1[0] < ga2[0]:
					ancPair = (ga1, ga2)
					modPair = (gm1, gm2)
				else:
					ancPair = ((ga2[0],-ga2[1]), (ga1[0],-ga1[1]))
					modPair = ((gm2[0],-gm2[1]), (gm1[0],-gm1[1]))
				subdicAncMod[ancPair].append((esp,modPair))
				dicModAnc[modPair].append( (anc,ancPair) )
				ga1 = ga2
				gm1 = gm2

	del todo[esp]
	print >> sys.stderr, "OK"

for esp in listSpecies:
	extractPairsFromSpecies(esp)

# Now that all the genomes have been loaded, let's empty the cache and restore intern
name_hash = {}
utils.myGenomes.intern = intern

print >> sys.stderr, "time for task1", time.time() - start
start = time.time()


# Returns all ancestral pairs designated by modern pairs
#########################################################
def getTargets(listAnc, lmodPair):
	lanc = collections.defaultdict(list)
	for modPair in lmodPair:
		t = dict(dicModAnc[modPair[1]])
		for otherAncName in listAnc:
			if otherAncName in t:
				lanc[(otherAncName, t[otherAncName])].append(modPair)
	return lanc.iteritems()
	

details = collections.defaultdict(lambda: collections.defaultdict(set))

def intersectAndPropagatePairs(anc):
	
	subdicAncMod = dicAncMod.pop(anc)
	# All pairs of child species, grouped by subtree
	pairs = [(x,subdicAncMod[x]) for (x,_) in phylTree.items[anc]]
	listAnc = phylTree.allDescendants[anc].intersection(phylTree.listAncestr)
	
	print >> sys.stderr, "Number of pairs for", anc, [(x[0],len(x[1])) for x in pairs]
	
	nbcons = 0
	# Iterate over pairs of subtrees
	for ((e1,gr1),(e2,gr2)) in itertools.combinations(pairs, 2):
		print >> sys.stderr, "Intersection between", e1, "and", e2, "...",
		# Iterate over the ancestral pairs found in both subtrees
		for ancPair in set(gr1).intersection(set(gr2)):
			nbcons += 1
			
			lmodPair1 = gr1[ancPair]    # Related modern pairs on the first subtree
			lmodPair2 = gr2[ancPair]    # Related modern pairs on the second subtree

			# And their reverse
			rlmodPair1 = [(x[0],revPair(x[1])) for x in lmodPair1]
			rlmodPair2 = [(x[0],revPair(x[1])) for x in lmodPair2]

			# Iterate over all the ancestors of the first subtree (incl. anc), and the modern pairs below them
			for ((ancName, ancPair), lmodPair) in getTargets(listAnc, lmodPair1):
				# Add evidence onto the ancestral pair of this ancestor:
				# - the subset of modern pairs of the first subtree
				# - all the modern pairs of the second subtree
				details[ancName][ancPair].update(lmodPair, lmodPair2)

			# Do the same for all the other subtree
			for ((ancName, ancPair), lmodPair) in getTargets(listAnc, lmodPair2):
				details[ancName][ancPair].update(lmodPair, lmodPair1)
	
			# And for the reverse pairs
			for ((ancName, ancPair), lmodPair) in getTargets(listAnc, rlmodPair1):
				details[ancName][ancPair].update(lmodPair, rlmodPair2)
	
			for ((ancName, ancPair), lmodPair) in getTargets(listAnc, rlmodPair2):
				details[ancName][ancPair].update(lmodPair, rlmodPair1)

		print >> sys.stderr, "OK"

	print >> sys.stderr, nbcons, "conserved pairs between descendants", anc

for anc in list(dicAncMod):
	intersectAndPropagatePairs(anc)

del dicModAnc
print >> sys.stderr, "time for task 2", time.time() - start
start = time.time()

# Results files.
def reportPairs(anc):
	pairs = details.pop(anc)

	# Accessory ancestor (required to compare against outgroups)
	if anc not in listAncestors:
		print >> sys.stderr, "Skipping", anc, "(not a target)"
		return

	print >> sys.stderr, len(pairs), "conserved pairs for", anc

	# -1 is the outgroup species, 1,2,3... are the descendantsdescendants
	ind = dict.fromkeys(set(phylTree.outgroupSpecies[anc]).intersection(listSpecies), -1)
	for (i,(x,_)) in enumerate(phylTree.items[anc]):
		ind.update(dict.fromkeys(set(phylTree.species[x]).intersection(listSpecies), i+1))
	
	res = arguments["OUT.pairwise"] % phylTree.fileName[anc]
	f = utils.myFile.openFile(res, "w")
	for ancPair in pairs:

		#assert ancPair[0][0] < ancPair[1][0]
		#assert set(pairs[ancPair]).isdisjoint([(x[0],revPair(x[1])) for x in pairs[ancPair]])

		# Calcul des poids
		weights = collections.defaultdict(int)
		for modPair in pairs[ancPair]:
			weights[ind[modPair[0]]] += 1
		weight = sum(x*y for (x,y) in itertools.combinations(weights.values(), 2))
	    
		print >> f, utils.myFile.myTSV.printLine(
			list(ancPair[0] + ancPair[1]) +
			[weight]
		)
	f.close()
		
for anc in list(details):
	reportPairs(anc)

print >> sys.stderr, "Elapsed time task3:", (time.time() - start), (time.time() - st)
