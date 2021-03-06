#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Convert a genome written as blocks of contigs to blocks of ancGenes
"""

import multiprocessing
import sys
import time

from joblib import Parallel, delayed

import utils.myFile
import utils.myGenomes
import utils.myGraph
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("target", str)],
    [("nbThreads", int, 0), ("IN.scaffoldsFile", str, ""), ("IN.contigsFile", str, ""), ("OUT.ancBlocksFile", str, "")],
    __doc__
)

# Load species tree - target ancestral genome and the extant species used to assemble blocs
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

def do(anc):
    (diags,singletons) = utils.myGraph.loadIntegr(arguments["IN.scaffoldsFile"] % phylTree.fileName[anc])

    print >> sys.stderr, "Loading reference contigs set from", arguments["IN.contigsFile"] % phylTree.fileName[anc], "...",
    ref = {}
    fi = utils.myFile.openFile(arguments["IN.contigsFile"] % phylTree.fileName[anc], "r")
    for (i,l) in enumerate(fi):
        ref[i+1] = l
    fi.close()
    print >> sys.stderr, "OK"

    print >> sys.stderr, "Writing ancBlocks of", anc, "...",
    lengths = []
    fo = utils.myFile.openFile(arguments["OUT.ancBlocksFile"] % phylTree.fileName[anc], "w")
    for (chrom,weights) in diags:
        li = []
        ls = []
        lw = []
        n = 0
        for (i,(c,s)) in enumerate(chrom):
            t = ref.pop(c)[:-1].split("\t")
            if i >= 1:
                lw.append(weights[i-1])
            n += len(t[2].split())
            if s > 0:
                li.append(t[2])
                ls.append(t[3])
                lw.append(t[4])
            else:
                li.extend(reversed(t[2].split()))
                ls.extend(-int(x) for x in reversed(t[3].split()))
                lw.extend(reversed(t[4].split()))

        lengths.append(n)
        print >> fo, utils.myFile.myTSV.printLine([t[0], n, utils.myFile.myTSV.printLine(li, delim=" "), utils.myFile.myTSV.printLine(ls, delim=" "), utils.myFile.myTSV.printLine(lw, delim=" ")])

    ns = 0
    for c in singletons:
        l = ref.pop(c)
        n = len(l[:-1].split("\t")[2].split())
        if n >= 2:
            lengths.append(n)
        else:
            ns += 1
        print >> fo, l,

    # S'assure que tous les contigs ont ete employes
    assert len(ref) == 0
    print >> sys.stderr, utils.myMaths.myStats.txtSummary(sorted(lengths)), "+", ns, "singletons"

start = time.time()
n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)
print >> sys.stderr, "Time elapsed:", time.time() - start
