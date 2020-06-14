#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Compute conserved block adjacencies across the given species set, and reconstruct blocks
    of blocks

    Usage:
        src/buildSynteny.integr-scaffolds.py example/data/Species.nwk A0 \
                -IN.ancBlocks=example/results/ancBlocks/denovo-all/blocks.%s.list.bz2 \
                -ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2 \
                -genesFiles=example/data/genes/genes.%s.list.bz2 \
                -OUT.ancBlocks=example/results/ancBlocks/denovo-all.scaffolds/blocks.%s.list.bz2 \
                -LOG.ancGraph=example/results/ancBlocks/denovo-all.scaffolds/graph.%s.txt.bz2
"""

import collections
import itertools
import multiprocessing
import sys
import time

import utils.myFile
import utils.myGenomes
import utils.myGraph
import utils.myMaths
import utils.myPhylTree
import utils.myTools
from utils.myTools import file

# Arguments
arguments = utils.myTools.checkArgs( \
    [("speciesTree", file), ("target", str)], \
    [("minimalWeight", int, 1), ("anchorSize", int, 2), ("minChromLength", int, 2), \
     ("nbThreads", int, 0),
     ("extantSpeciesFilter", str, ""), \
     ("IN.ancBlocks", str, ""), \
     ("LOG.ancGraph", str, ""),
     ("OUT.ancBlocks", str, ""), \
     ("genesFiles", str, ""), \
     ("ancGenesFiles", str, "")], \
    __doc__ \
    )


def rewriteGenome(genome, dic):
    newGenome = {}
    for chrom in genome.chrList[utils.myGenomes.ContigType.Chromosome] + genome.chrList[
        utils.myGenomes.ContigType.Scaffold]:
        tmp = [(dic[(chrom, i)], gene.strand) for (i, gene) in enumerate(genome.lstGenes[chrom]) if (chrom, i) in dic]
        tmp = [(i, s * sa) for ((i, sa), s) in tmp]
        if len(tmp) > 0:
            newGenome[chrom] = [i for (i, l) in itertools.groupby(tmp)]
    return newGenome


def getExtremities(genome):
    extr1 = {}
    extr2 = {}
    for (chrom, l) in genome.items():
        (i0, s0) = l[0]
        (i1, s1) = l[-1]
        extr1[(i0, s0)] = (chrom, 1)
        extr2[(i1, s1)] = (chrom, 1)
        extr1[(i1, -s1)] = (chrom, -1)
        extr2[(i0, -s0)] = (chrom, -1)
    return (extr1, extr2)


def getAllAdj(anc, dicGenomesAnc):
    allAdj = {}
    for esp in listSpecies:

        dicA = {}
        dicM = {}
        stats = []
        print("Gene order comparison between %s and %s ..." % (anc, esp), end=' ', file=sys.stderr)

        for (n, ((c1, d1), (c2, d2), da)) in enumerate(
                utils.myGraph.calcDiags(dicGenomes[esp], dicGenomesAnc, genesAnc[phylTree.dicParents[anc][esp]],
                                        orthosFilter=utils.myGraph.OrthosFilterType.InBothSpecies,
                                        minChromLength=arguments["minChromLength"])):
            if len(da) < arguments["anchorSize"]:
                continue
            # print "DIAG", anc, esp, n, (c1,c2), len(da), (d1,d2,da)
            for ((i1, s1), (i2, s2)) in zip(d1, d2):
                dicM[(c1, i1)] = (n, s1)
                dicA[(c2, i2)] = (n, s1)
            stats.append(len(da))

        print(utils.myMaths.myStats.txtSummary(stats), end=' ', file=sys.stderr)

        newGA = rewriteGenome(dicGenomesAnc, dicA)
        # List of selected blocks in the ancestor
        notdup = set()
        for cA in newGA:
            notdup.update(x[0] for x in newGA[cA])
        # print anc, esp, "ANC", cA, len(newGA[cA]), newGA[cA]
        newGM = rewriteGenome(dicGenomes[esp], dicM)

        (extr1, extr2) = getExtremities(newGA)

        tmp = []
        for (cM, l) in newGM.items():
            # print anc, esp, "MOD", cM, len(newGM[cM]), newGM[cM]
            # Allows to select during a segmental duplication, the same block than the ancestor
            l = [x for x in l if x[0] in notdup]
            # print anc, esp, "FMOD", cM, len(l), l
            for (x1, x2) in utils.myTools.myIterator.slidingTuple(l):
                if (x1 in extr2) and (x2 in extr1):
                    tmp.append((extr2[x1], extr1[x2]))
                # print "ADJ", anc, esp, (x1,x2), (extr2[x1],extr1[x2])

        allAdj[esp] = set(tmp)
        allAdj[esp].update(((i2, -s2), (i1, -s1)) for ((i1, s1), (i2, s2)) in tmp)
        print("%d adjacencies / %d blocks" % (len(tmp), len(newGA)), file=sys.stderr)
    return allAdj


def do(anc):
    # Redirect the standard output to a file
    ini_stdout = sys.stdout
    sys.stdout = utils.myFile.openFile(arguments["LOG.ancGraph"] % phylTree.fileName[anc], "w")

    dicGenomesAnc = utils.myGenomes.Genome(arguments["IN.ancBlocks"] % phylTree.fileName[anc], ancGenes=genesAnc[anc],
                                             withDict=False)

    allAdj = getAllAdj(anc, dicGenomesAnc)

    gr = utils.myGraph.WeightedDiagGraph()
    for (e1, e2) in toStudy[anc]:
        for x in allAdj[e1] & allAdj[e2]:
            gr.addDiag(x)
    # gr.printIniGraph()
    del allAdj
    gr.cleanGraphTopDown(2 * arguments["minimalWeight"])

    stats = []
    f = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    notseen = set(dicGenomesAnc.lstGenes)

    def toString(x, rev=False):
        lg = [genesAnc[anc].dicGenes[gene.names[0]].index for gene in dicGenomesAnc.lstGenes[x]]
        ls = [gene.strand for gene in dicGenomesAnc.lstGenes[x]]
        if rev:
            lg.reverse()
            ls = [-x for x in reversed(ls)]
        return (utils.myFile.myTSV.printLine(lg, delim=" "), utils.myFile.myTSV.printLine(ls, delim=" "))

    # First, all the fusions of blocks are written
    for (res, scores) in gr.getBestDiags():
        if len(res) == 1:
            continue
        l = 0
        tg = []
        ts = []
        tw = []
        for (x, s) in res:
            notseen.remove(x)
            tmp = toString(x, rev=(s < 0))
            tg.append(tmp[0])
            ts.append(tmp[1])
            if s > 0:
                tw.extend(dicGenomesAnc.support[x])
            else:
                tw.extend(reversed(dicGenomesAnc.support[x]))
            if len(scores) > 0:
                tw.append(str(scores.pop(0)))
            l += len(dicGenomesAnc.lstGenes[x])
        print(utils.myFile.myTSV.printLine([anc, l, " ".join(tg), " ".join(ts), " ".join(tw)]), file=f)
        stats.append(l)

    # Write the single blocks
    for x in sorted(notseen):
        print(utils.myFile.myTSV.printLine(
            (anc, len(dicGenomesAnc.lstGenes[x])) + toString(x) + ("(%d)" % len(dicGenomesAnc.lstGenes[x]),)), file=f)
        if len(dicGenomesAnc.lstGenes[x]) > 1:
            stats.append(len(dicGenomesAnc.lstGenes[x]))

    print("Output blocks of", anc, utils.myMaths.myStats.txtSummary(stats), "+", len(
        genesAnc[anc].lstGenes[None]) - sum(stats), "singletons", file=sys.stderr)
    f.close()

    # Revert to the true standard output
    sys.stdout.close()
    sys.stdout = ini_stdout


# Load species tree - target ancestral genome and the extant species used to assemble blocks
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

(listSpecies, targets, accessoryAncestors) = phylTree.getTargetsForPairwise(arguments["target"], arguments["extantSpeciesFilter"])
listSpecies = sorted(listSpecies)

# Gene names don't matter by themselves. What is important is that they link the dicGenomes
# and the genesAnc. utils.myGenomes uses sys.intern() to make them refer to the same string
# instances to reduce the memory usage. Here we go further and we replace the names with
# unique integers (the smallest object in Python) to save even more memory
name_hash = {}
n_names = 0
def myintern(s):
    t = name_hash.get(s)
    if t:
        return t
    global n_names
    n_names += 1
    name_hash[s] = n_names
    return n_names

utils.myGenomes.intern = myintern

genesAnc = {}
for anc in sorted(targets.union(accessoryAncestors)):
    genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])

# Here's another trick. Once we have loaded all the ancestral genes, we have indexed the
# names of all extant genes. Since extant genes are unique, we can clear the names used by
# each species right after having loaded it
# Note: this set will also contain the chromosome names
species_names = set()
def myintern_species(s):
    species_names.add(s)
    return myintern(s)

utils.myGenomes.intern = myintern_species

dicGenomes = {}
for e in listSpecies:
    dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[e], withDict=False)
    for s in species_names:
        del name_hash[s]
    species_names = set()

# Now that the names have been replaced, let's empty the cache and restore intern
name_hash = {}
utils.myGenomes.intern = sys.intern

toStudy = collections.defaultdict(list)
for (e1, e2) in itertools.combinations(listSpecies, 2):
    for anc in targets.intersection(phylTree.dicLinks[e1][e2][1:-1]):
        toStudy[anc].append((e1, e2))

start = time.time()
n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
multiprocessing.Pool(n_cpu).map(do, sorted(targets))
print("Elapsed time:", (time.time() - start), file=sys.stderr)
