#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Build the ancestral graph of conserved gene adjacencies (inferred from pairwise comparisons) and linearize it.

    Usage:
        src/buildSynteny.integr-denovo.py example/data/Species.nwk A0 \
                example/results/pairwise/pairs-all/%s.list.bz2 \
                -ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2 \
                -OUT.ancBlocks=example/results/ancBlocks/denovo-all/blocks.%s.list.bz2 \
                -LOG.ancGraph=example/results/ancBlocks/denovo-all/graph.%s.txt.bz2
"""

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
arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("target", str), ("pairwise", str)],
    [("minimalWeight", int, 1), ("searchLoops", bool, False),
     ("LOG.ancGraph", str, ""),
     ("OUT.ancBlocks", str, ""),
     ("ancGenesFiles", str, ""),
     ("nbThreads", int, 0),
     ],
    __doc__
)


def do(anc):

    diags = utils.myGraph.loadConservedPairsAnc(arguments["pairwise"] % phylTree.fileName[anc])

    g = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc], withDict=False).lstGenes
    singletons = set(range(len(g[None]))) if None in g else set(g)

    print("Blocks of %s ..." % anc, end=' ', file=sys.stderr)

    graph = utils.myGraph.WeightedDiagGraph()
    for x in diags:
        graph.addLink(*x)

    # Redirect the standard output to a file
    ini_stdout = sys.stdout
    sys.stdout = utils.myFile.openFile(arguments["LOG.ancGraph"] % phylTree.fileName[anc], "w")

    graph.printIniGraph()

    # Cut the graph in subgraph
    graph.cleanGraphTopDown(arguments["minimalWeight"], searchLoops=arguments["searchLoops"])

    f = utils.myFile.openFile(arguments["OUT.ancBlocks"] % phylTree.fileName[anc], "w")
    s = []

    # Graph linearisation
    for (d, dw) in graph.getBestDiags():

        if len(d) == 1:
            continue

        ds = [x[1] for x in d]
        da = [x[0] for x in d]

        s.append(len(da))
        singletons.difference_update(da)
        res = [anc, len(da), utils.myFile.myTSV.printLine(da, " "), utils.myFile.myTSV.printLine(ds, " "),
               utils.myFile.myTSV.printLine(dw, " ")]
        print(utils.myFile.myTSV.printLine(res), file=f)

    for x in singletons:
        print(utils.myFile.myTSV.printLine([anc, 1, x, 1, ""]), file=f)
    f.close()
    print("OK", file=sys.stderr)
    print(anc,  utils.myMaths.myStats.syntheticTxtSummary(s), "+ %d singletons OK" % len(singletons), file=sys.stderr)

    # Revert to the true standard output
    sys.stdout.close()
    sys.stdout = ini_stdout


start = time.time()

# Load species tree and target ancestral genome
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

print("Targets:", sorted(targets), file=sys.stderr)

n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
multiprocessing.Pool(n_cpu).map(do, sorted(targets))

print("Elapsed time:", (time.time() - start), file=sys.stderr)
