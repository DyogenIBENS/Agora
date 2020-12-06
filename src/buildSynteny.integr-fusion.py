#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v1.5
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Rallonge et fusionne des blocs integres de base, grace a des blocs de syntenie pairwise relaxes
"""

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
    [("phylTree.conf", file), ("target", str), ("pairwiseDiags", str)],
    [("minimalWeight", int, 1), ("searchLoops", bool, True), ("onlySingletons", bool, False),
     ("nbThreads", int, 0),
     ("IN.ancBlocks", str, ""), ("OUT.ancBlocks", str, ""), ("LOG.ancGraph", str, "extend_log/%s.log.bz2")],
    __doc__
)


def do(anc):
    # Redirect the standard output to a file
    ini_stdout = sys.stdout
    sys.stdout = utils.myFile.openFile(arguments["LOG.ancGraph"] % phylTree.fileName[anc], "w")

    graph = utils.myGraph.WeightedDiagGraph()
    (integr, singletons) = utils.myGraph.loadIntegr(arguments["IN.ancBlocks"] % phylTree.fileName[anc])
    if not arguments["onlySingletons"]:
        for (b, w) in integr:
            graph.addWeightedDiag(b, [x + 10000 for x in w])

    diags = utils.myGraph.loadConservedPairsAnc(arguments["pairwiseDiags"] % phylTree.fileName[anc])
    for x in diags:
        if not arguments["onlySingletons"] or ((x[0][0] in singletons) and (x[1][0] in singletons)):
            graph.addLink(*x)

    print >> sys.stderr, "Blocs integres de %s ..." % anc,

    # cutting the graph
    graph.cleanGraphTopDown(arguments["minimalWeight"], searchLoops=arguments["searchLoops"])

    f = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    s = []

    if arguments["onlySingletons"]:
        # If this option is set, the blocks are printed as they are
        for (b, w) in integr:
            print >> f, utils.myFile.myTSV.printLine([anc, len(b),
                                                      utils.myFile.myTSV.printLine([x[0] for x in b], " "),
                                                      utils.myFile.myTSV.printLine([x[1] for x in b], " "),
                                                      utils.myFile.myTSV.printLine(w, " ")])

    # Integrated blocks extraction
    for (d, dw) in graph.getBestDiags():

        if len(d) == 1:
            continue

        da = [x[0] for x in d]
        ds = [x[1] for x in d]

        s.append(len(da))
        singletons.difference_update(da)
        print >> f, utils.myFile.myTSV.printLine(
            [anc, len(da), utils.myFile.myTSV.printLine(da, " "), utils.myFile.myTSV.printLine(ds, " "),
             utils.myFile.myTSV.printLine(dw, " ")])

    for x in singletons:
        print >> f, utils.myFile.myTSV.printLine([anc, 1, x, 1, ""])
    f.close()
    print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "+ %d singletons OK" % len(singletons)

    # Revert to the true standard output
    sys.stdout.close()
    sys.stdout = ini_stdout


start = time.time()


# Load species tree and target ancestral genome
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
targets = phylTree.getTargetsAnc(arguments["target"])

n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)

print >> sys.stderr, "Elapsed time:", (time.time() - start)
