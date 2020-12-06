#!/usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v1.5
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Transform reconstructed ancestral lists of blocks in formated ancestral genomes (lists of genes)

    usage:
        src/convert.ancGenomes.blocks-to-genes.py example/data/Species.conf A0 \
                -IN.ancBlocks=example/results/integrDiags/final/diags.%s.list.bz2 \
                -OUT.ancGenomes=example/results/ancGenomes/final/ancGenome.%s.list.bz2 \
                -ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2
"""

import multiprocessing
import sys
import time

from joblib import Parallel, delayed

import utils.myGenomes
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("target", str)],
    [("nbThreads", int, 0), ("IN.ancBlocks", str, ""), ("ancGenesFiles", str, ""), ("OUT.ancGenomes", str, "ancGenomes/ancGenome.%s.list.bz2")],
    __doc__
)

# Load species tree - target ancestral genome and the extant species used to assemble blocs
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

def do(anc):
    ancGenes = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])
    genome = utils.myGenomes.Genome(arguments["IN.ancBlocks"] % phylTree.fileName[anc], ancGenes=ancGenes)
    ancGenomeFile = utils.myFile.openFile(arguments["OUT.ancGenomes"] % phylTree.fileName[anc], "w")
    genome.printEnsembl(ancGenomeFile)
    ancGenomeFile.close()

start = time.time()
n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)
print >> sys.stderr, "Time elapsed:", time.time() - start
