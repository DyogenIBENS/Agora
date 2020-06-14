#!/usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Convert an ancestral genome written (blocks of genes) to the BED format

    Usage:
        src/convert.ancGenomes.blocks-to-genes.py example/data/Species.nwk A0 \
                -IN.ancBlocks=example/results/best-pass1.best-pass2/blocks.%s.list.bz2 \
                -OUT.ancGenomes=example/results/ancGenomes/generic-workflow/ancGenome.%s.list.bz2 \
                -ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2
"""

import itertools
import multiprocessing
import sys
import time

import utils.myFile
import utils.myGenomes
import utils.myPhylTree
import utils.myTools
from utils.myTools import file

# Arguments
arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("target", str)],
    [("IN.ancBlocks", str, ""), ("ancGenesFiles", str, ""), ("OUT.ancGenomes", str, ""),
     ("nbThreads", int, 0), ("orderBySize", bool, False),
    ],
    __doc__
)

# Load species tree - target ancestral genome and the extant species used to assemble blocks
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

def do(anc):
    genome = utils.myGenomes.Genome(arguments["IN.ancBlocks"] % phylTree.fileName[anc], withDict=False)

    block_names = list(genome.lstGenes.keys())
    names = {}
    if arguments["orderBySize"]:
        block_names.sort(key=lambda c: len(genome.lstGenes[c]), reverse=True)
        n_CARs = 0
        n_singletons = 0
        for s in block_names:
            if len(genome.lstGenes[s]) == 1:
                n_singletons += 1
                names[s] = "singleton_%d" % n_singletons
            else:
                n_CARs += 1
                names[s] = "CAR_%d" % n_CARs
    else:
        # Simply use integers
        for (i, s) in enumerate(block_names):
            names[s] = i

    # Print the ancestral genome, using the actual gene names provided by ancGenesFiles
    ancGenes = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc]).lstGenes[None]
    ancGenomeFile = utils.myFile.openFile(arguments["OUT.ancGenomes"] % phylTree.fileName[anc], "w")
    for s in block_names:
        for gene in genome.lstGenes[s]:
            print(utils.myFile.myTSV.printLine([names[s], gene.beginning, gene.end, gene.strand, " ".join(ancGenes[gene.names[0]].names)]), file=ancGenomeFile)
    ancGenomeFile.close()

start = time.time()
n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
multiprocessing.Pool(n_cpu).map(do, sorted(targets))
print("Time elapsed:", time.time() - start, file=sys.stderr)
