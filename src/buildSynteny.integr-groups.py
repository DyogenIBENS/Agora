#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Combine previous blocks when they are side by side in extant species.
    This could be compare to an "assembly" step of previous blocs.

        Usage:
            src/buildSynteny.integr-groups.py Species.conf Boreoeutheria Boreoeutheria -IN.ancDiags=diags/integr/denovo-all/anc/diags.%s.list.bz2 -OUT.ancDiags=diags/integr/denovo-all.groups/anc/diags.%s.list.bz2
            -ancGenesFiles=ancGenes/all/ancGenes.%s.list.bz2 -genesFiles=genes.%s.list.bz2 2> diags/integr/denovo-all.groups/log
"""

import collections
import itertools
import sys

import utils.myFile
import utils.myGenomes
import utils.myGraph
import utils.myMaths
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs( \
    [("phylTree.conf", file), ("target", str), ("usedSpecies", str)], \
    [("minimalWeight", int, 1), ("anchorSize", int, 2), ("minChromLength", int, 2), \
     ("IN.ancDiags", str, ""), \
     ("LOG.ancGraph", str, "groups_log/%s.log.bz2"),
     ("OUT.ancDiags", str, "anc/diags.%s.list.bz2"), \
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
    for (chrom, l) in genome.iteritems():
        (i0, s0) = l[0]
        (i1, s1) = l[-1]
        extr1[(i0, s0)] = (chrom, 1)
        extr2[(i1, s1)] = (chrom, 1)
        extr1[(i1, -s1)] = (chrom, -1)
        extr2[(i0, -s0)] = (chrom, -1)
    return (extr1, extr2)


def getAllAdj(anc):
    allAdj = {}
    for esp in listSpecies:

        dicA = {}
        dicM = {}
        stats = []
        print >> sys.stderr, "Diagonals extraction between %s and %s ..." % (anc, esp),

        for (n, ((c1, d1), (c2, d2), da)) in enumerate(
                utils.myGraph.calcDiags(dicGenomes[esp], dicGenomes[anc], genesAnc[phylTree.dicParents[anc][esp]],
                                        orthosFilter=utils.myGraph.OrthosFilterType.InBothSpecies,
                                        minChromLength=arguments["minChromLength"])):
            if len(da) < arguments["anchorSize"]:
                continue
            # print "DIAG", anc, esp, n, (c1,c2), len(da), (d1,d2,da)
            for ((i1, s1), (i2, s2)) in itertools.izip(d1, d2):
                dicM[(c1, i1)] = (n, s1)
                dicA[(c2, i2)] = (n, s1)
            stats.append(len(da))

        print >> sys.stderr, utils.myMaths.myStats.txtSummary(stats),

        newGA = rewriteGenome(dicGenomes[anc], dicA)
        # List of selected blocks in the ancestor
        notdup = set()
        for cA in newGA:
            notdup.update(x[0] for x in newGA[cA])
        # print anc, esp, "ANC", cA, len(newGA[cA]), newGA[cA]
        newGM = rewriteGenome(dicGenomes[esp], dicM)

        (extr1, extr2) = getExtremities(newGA)

        tmp = []
        for (cM, l) in newGM.iteritems():
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
        print >> sys.stderr, "%d adjacencies / %d blocs" % (len(tmp), len(newGA))
    return allAdj


def do(anc):
    allAdj = getAllAdj(anc)

    gr = utils.myGraph.WeightedDiagGraph()
    for (e1, e2) in toStudy[anc]:
        for x in allAdj[e1] & allAdj[e2]:
            gr.addDiag(x)
    # gr.printIniGraph()
    gr.cleanGraphTopDown(2 * arguments["minimalWeight"])

    stats = []
    f = utils.myFile.openFile(arguments["OUT.ancDiags"] % phylTree.fileName[anc], "w")
    notseen = set(dicGenomes[anc].lstGenes)

    def toString(x, rev=False):
        lg = [genesAnc[anc].dicGenes[gene.names[0]].index for gene in dicGenomes[anc].lstGenes[x]]
        ls = [gene.strand for gene in dicGenomes[anc].lstGenes[x]]
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
            tw.append("(%d)" % len(dicGenomes[anc].lstGenes[x]))
            if len(scores) > 0:
                tw.append(str(scores.pop(0)))
            l += len(dicGenomes[anc].lstGenes[x])
        print >> f, utils.myFile.myTSV.printLine([anc, l, " ".join(tg), " ".join(ts), " ".join(tw)])
        stats.append(l)

    # Write the single blocks
    for x in sorted(notseen):
        print >> f, utils.myFile.myTSV.printLine(
            (anc, len(dicGenomes[anc].lstGenes[x])) + toString(x) + ("(%d)" % len(dicGenomes[anc].lstGenes[x]),))
        if len(dicGenomes[anc].lstGenes[x]) > 1:
            stats.append(len(dicGenomes[anc].lstGenes[x]))

    print >> sys.stderr, "Integrated blocs of", anc, utils.myMaths.myStats.txtSummary(stats), "+", len(
        genesAnc[anc].lstGenes[None]) - sum(stats), "singletons"
    f.close()


# Load species tree - target ancestral genome and the extant species used to assemble blocs
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

targets = phylTree.getTargetsAnc(arguments["target"])
listSpecies = phylTree.getTargetsSpec(arguments["usedSpecies"])

dicGenomes = {}
for e in listSpecies:
    dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFiles"] % phylTree.fileName[e], withDict=False)

genesAnc = {}
for anc in targets:
    genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])
    dicGenomes[anc] = utils.myGenomes.Genome(arguments["IN.ancDiags"] % phylTree.fileName[anc], ancGenes=genesAnc[anc],
                                             withDict=False)

for anc in [phylTree.dicParents[e][a] for (e, a) in itertools.product(listSpecies, targets)]:
    if anc not in genesAnc:
        genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])

toStudy = collections.defaultdict(list)
for (e1, e2) in itertools.combinations(listSpecies, 2):
    for anc in targets.intersection(phylTree.dicLinks[e1][e2][1:-1]):
        toStudy[anc].append((e1, e2))

for anc in targets:
    # Redirect the standard output to a file
    ini_stdout = sys.stdout
    sys.stdout = utils.myFile.openFile(arguments["LOG.ancGraph"] % phylTree.fileName[anc], "w")
    do(anc)
    # Revert to the true standard output
    sys.stdout.close()
    sys.stdout = ini_stdout
