#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Builds the ancestral graph of gene adjacencies infered from pairwise comparisons and linearizes it.

    usage:
        src/buildSynteny.integr-denovo.py Species.conf Boreoeutheria diags/pairwise/pairs-all/%s.list.bz2
        -ancGenesFiles=ancGenes/all/ancGenes.%s.list.bz2 -OUT.ancDiags=diags/integr/denovo-all/anc/diags.%s.list.bz2
        2> diags/integr/denovo-all/log
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
    [("minimalWeight", int, 1), ("searchLoops", bool, False),
     ("LOG.ancGraph", str, "denovo_log/%s.log.bz2"),
     ("OUT.ancDiags", str, "anc/diags.%s.list.bz2"),
     ("ancGenesFiles", str, ""),
     ("nbThreads", int, 0),
     ],
    __doc__
)


def do(anc):

    diags = utils.myGraph.loadConservedPairsAnc(arguments["pairwiseDiags"] % phylTree.fileName[anc])

    g = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc], withDict=False).lstGenes
    singletons = set(xrange(len(g[None]))) if None in g else set(g)

    print >> sys.stderr, "Integrated blocs of %s ..." % anc,

    graph = utils.myGraph.WeightedDiagGraph()
    for x in diags:
        graph.addLink(*x)

    # Redirect the standard output to a file
    ini_stdout = sys.stdout
    sys.stdout = utils.myFile.openFile(arguments["LOG.ancGraph"] % phylTree.fileName[anc], "w")

    graph.printIniGraph()

    # Cut the graph in subgraph
    graph.cleanGraphTopDown(arguments["minimalWeight"], searchLoops=arguments["searchLoops"])

    f = utils.myFile.openFile(arguments["OUT.ancDiags"] % phylTree.fileName[anc], "w")
    s = []

    # Graph linearizing
    for (d, dw) in graph.getBestDiags():

        if len(d) == 1:
            continue

        ds = [x[1] for x in d]
        da = [x[0] for x in d]

        s.append(len(da))
        singletons.difference_update(da)
        res = [anc, len(da), utils.myFile.myTSV.printLine(da, " "), utils.myFile.myTSV.printLine(ds, " "),
               utils.myFile.myTSV.printLine(dw, " ")]
        print >> f, utils.myFile.myTSV.printLine(res)

    for x in singletons:
        print >> f, utils.myFile.myTSV.printLine([anc, 1, x, 1, ""])
    f.close()
    print >> sys.stderr, "OK"
    print >> sys.stderr, anc,  utils.myMaths.myStats.syntheticTxtSummary(s), "+ %d singletons OK" % len(singletons)

    # Revert to the true standard output
    sys.stdout.close()
    sys.stdout = ini_stdout


start = time.time()

# Load species tree and target ancestral genome
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
targets = phylTree.getTargetsAnc(arguments["target"])

print >> sys.stderr, targets

n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)

print >> sys.stderr, "Elapsed time:", (time.time() - start)
