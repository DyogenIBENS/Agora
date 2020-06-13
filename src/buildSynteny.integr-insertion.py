#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Insere les blocs de genes non contraints a l'interieur et aux extremites
"""

import collections
import itertools
import multiprocessing
import sys
import time

from joblib import Parallel, delayed

import utils.myFile
import utils.myGenomes
import utils.myGraph
import utils.myMaths
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("target", str), ("pairwise", str)],
    [("IN.ancBlocks", str, ""), ("OUT.ancBlocks", str, ""), ("REF.ancBlocks", str, ""), ("LOG.ancGraph", str, "halfinsert_log/%s.log.bz2"),
     ("nbThreads", int, 0),
     ("selectionFunction", str, "newscore/float(oldscore) if oldscore else newscore")],
    __doc__
)


# reverse gene
###############
def rev(xxx_todo_changeme):
    (g, s) = xxx_todo_changeme
    return (g, -s)


selectionFunction = eval("lambda newscore, oldscore: " + arguments["selectionFunction"])


def do(anc):
    # Redirect the standard output to a file
    ini_stdout = sys.stdout
    sys.stdout = utils.myFile.openFile(arguments["LOG.ancGraph"] % phylTree.fileName[anc], "w")

    pairwiseDiags = loadPairwise(anc)
    (integr, singletons) = utils.myGraph.loadIntegr(arguments["IN.ancBlocks"] % phylTree.fileName[anc])
    (_, refsing) = utils.myGraph.loadIntegr(arguments["REF.ancBlocks"] % phylTree.fileName[anc])

    for (x, _) in integr:
        assert (x[0][0] in refsing) == (x[-1][0] in refsing), x
    iniblocks = [x for x in integr if x[0][0][0] not in refsing]
    toaddblocks = [x for x in integr if x[0][0][0] in refsing]
    assert singletons.issubset(refsing)
    for x in sorted(singletons):
        toaddblocks.append(([(x, 1)], []))
    print("iniblocks", len(iniblocks))
    print("toaddblocks", len(toaddblocks))

    extr = {}
    for (i, (lg, _)) in enumerate(toaddblocks):
        extr[lg[0]] = (i, True)
        extr[rev(lg[-1])] = (i, False)
    # print "extr", len(extr), extr.keys()[:10]

    # Possible links inside the genes intervals
    possins = []
    for (i, (lg, lw)) in enumerate(iniblocks):
        for (j, ((start, end), weight)) in enumerate(zip(utils.myTools.myIterator.slidingTuple(lg), lw)):

            def filt(items):
                return [x for x in items if (x[0] in extr) and (x[1] > weight)]

            out1 = filt(pairwiseDiags[start].items())
            out2 = filt(pairwiseDiags[rev(end)].items())

            for x in out1:
                (k, f) = extr[x[0]]
                possins.append((selectionFunction(x[1], weight), (i, j, "left"), k, not f))

            for x in out2:
                (k, f) = extr[x[0]]
                possins.append((selectionFunction(x[1], weight), (i, j, "right"), k, f))

    possins.sort(reverse=True)

    # Best links selection
    dicIntervBlock = {}
    dicBlockInterv = {}
    for x in possins:
        if (x[1] not in dicIntervBlock) and (x[2] not in dicBlockInterv):
            (i, j, _) = x[1]
            k = x[2]
            toreverse = x[3]
            if toreverse:
                # print "reverse block", k
                toaddblocks[k] = ([rev(g) for g in reversed(toaddblocks[k][0])], list(reversed(toaddblocks[k][1])))
            dicIntervBlock[x[1]] = x[2]
            dicBlockInterv[x[2]] = x[1]
            del extr[toaddblocks[x[2]][0][0]]
            del extr[rev(toaddblocks[x[2]][0][-1])]

    # Possible links at the end of a block
    possadd = []
    for (i, (lg, _)) in enumerate(iniblocks):
        possadd.extend((s, (i, "end"), extr[x]) for (x, s) in pairwiseDiags[lg[-1]].items() if x in extr)
        possadd.extend((s, (i, "start"), extr[x]) for (x, s) in pairwiseDiags[rev(lg[0])].items() if x in extr)
    possadd.sort(reverse=True)

    # Best link selection
    dicExtrBlock = {}
    dicBlockExtr = {}
    for x in possadd:
        if (x[1] not in dicExtrBlock) and (x[2][0] not in dicBlockExtr):
            k = x[2][0]
            if x[2][1] == (x[1][1] == "start"):
                # print "reverse block", k
                toaddblocks[k] = ([rev(g) for g in reversed(toaddblocks[k][0])], list(reversed(toaddblocks[k][1])))
            dicExtrBlock[x[1]] = x[2][0]
            dicBlockExtr[x[2][0]] = x[1]
            del extr[toaddblocks[x[2][0]][0][0]]
            del extr[rev(toaddblocks[x[2][0]][0][-1])]

    # Add a bloc at the end of the current one
    def addAfter(k, newb, neww):
        # print "add block", k, len(toaddblocks[k][0]), toaddblocks[k]
        # print "start", newb[-1], "->", toaddblocks[k][0][0], "|", pairwiseDiags[newb[-1]].get( toaddblocks[k][0][0] ,0)
        # print "end", toaddblocks[k][0][-1]
        assert len(toaddblocks[k][0]) >= 1
        neww.append(pairwiseDiags[newb[-1]].get(toaddblocks[k][0][0], 0))
        newb.extend(toaddblocks[k][0])
        neww.extend(toaddblocks[k][1])
        toaddblocks[k] = ([], [])

    #print >> sys.stderr, "Final blocks of ", anc,
    print("Final blocks of ", anc, end=' ')
    f = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    ll = []
    # Build and print new chromosomes
    for (i, (inib, iniw)) in enumerate(iniblocks):
        assert len(inib) >= 2

        k = dicExtrBlock.get((i, "start"))
        if k is None:
            newb = []
            neww = []
        else:
            assert len(toaddblocks[k][0]) >= 1
            newb = toaddblocks[k][0]
            neww = toaddblocks[k][1] + [pairwiseDiags[newb[-1]][inib[0]]]
            # print "chromosome start", k, len(newb), newb
            # print "end", newb[-1], neww[-1]
            toaddblocks[k] = ([], [])

        for (j, (g, weight)) in enumerate(zip(inib, iniw)):
            # print "initial interval", "%d/%d" % (i,j), g, "->", inib[j+1] if j+1 < len(inib) else None, "|", weight
            newb.append(g)

            kl = dicIntervBlock.get((i, j, "left"))
            kr = dicIntervBlock.get((i, j, "right"))
            if (kl is None) and (kr is None):
                neww.append(weight)
            else:
                if kl is not None:
                    # print "left"
                    addAfter(kl, newb, neww)
                if kr is not None:
                    # print "right"
                    addAfter(kr, newb, neww)
                if j + 1 < len(inib):
                    # print "added weight", pairwiseDiags[newb[-1]].get( inib[j+1] ,0)
                    neww.append(pairwiseDiags[newb[-1]].get(inib[j + 1], 0))

        newb.append(inib[-1])
        k = dicExtrBlock.get((i, "end"))
        if k is not None:
            # print "chromosome end"
            addAfter(k, newb, neww)

        assert len(newb) >= 2
        assert len(newb) == (len(neww) + 1)
        ll.append(len(newb))
        print(utils.myFile.myTSV.printLine([anc, len(newb),
                                                  utils.myFile.myTSV.printLine([x[0] for x in newb], delim=" "),
                                                  utils.myFile.myTSV.printLine([x[1] for x in newb], delim=" "),
                                                  utils.myFile.myTSV.printLine(neww, delim=" ")]), file=f)

    sing = 0
    for (newb, news) in toaddblocks:
        if len(newb) == 0:
            continue
        if len(newb) > 1:
            ll.append(len(newb))
        else:
            sing += len(newb)
        print(utils.myFile.myTSV.printLine([anc, len(newb),
                                                  utils.myFile.myTSV.printLine([x[0] for x in newb], delim=" "),
                                                  utils.myFile.myTSV.printLine([x[1] for x in newb], delim=" "),
                                                  utils.myFile.myTSV.printLine(news, delim=" ")]), file=f)
    f.close()
    print(utils.myMaths.myStats.txtSummary(ll), "+", sing, "singletons", file=sys.stderr)

    # Revert to the true standard output
    sys.stdout.close()
    sys.stdout = ini_stdout

start = time.time()
# Load species tree and target ancestral genome

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])


def loadPairwise(anc):
    lstPairwise = utils.myGraph.loadConservedPairsAnc(arguments["pairwise"] % phylTree.fileName[anc])
    pairwiseDiags = collections.defaultdict(dict)
    for d in lstPairwise:
        pairwiseDiags[d[0]][d[1]] = d[2]
        pairwiseDiags[rev(d[1])][rev(d[0])] = d[2]
    return pairwiseDiags


n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)

print("Elapsed time:", (time.time() - start), file=sys.stderr)
