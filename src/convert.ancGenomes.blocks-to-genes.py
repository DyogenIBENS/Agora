#!/usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Transform reconstructed ancestral lists of blocks in formated ancestral genomes (lists of genes)

    usage:
        src/convert.ancGenomes.blocks-to-genes.py example/data/Species.conf A0 \
                -IN.ancBlocks=example/results/ancBlocks/final/diags.%s.list.bz2 \
                -OUT.ancGenomes=example/results/ancGenomes/final/ancGenome.%s.list.bz2 \
                -ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2
"""

import itertools
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
    [("IN.ancBlocks", str, ""), ("IN.pairwise", str, ""), ("ancGenesFiles", str, ""), ("OUT.ancGenomes", str, "ancGenomes/ancGenome.%s.list.bz2"),
     ("nbThreads", int, 0), ("orderBySize", bool, False),
    ],
    __doc__
)

# Load species tree - target ancestral genome and the extant species used to assemble blocs
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
targets = phylTree.getTargetsAnc(arguments["target"])

def do(anc):
    genome = utils.myGenomes.Genome(arguments["IN.ancBlocks"] % phylTree.fileName[anc], withDict=False)

    block_names = genome.lstGenes.keys()
    if arguments["orderBySize"]:
        block_names.sort(key=lambda c: len(genome.lstGenes[c]), reverse=True)

    # If pairwise are available, load them and annotate the blocks as contigs and scaffolds
    names = {}
    if arguments["IN.pairwise"]:
        # Load the file and record the pairs in a set
        diags = utils.myGraph.loadConservedPairsAnc(arguments["IN.pairwise"] % phylTree.fileName[anc])
        pairwise = set()
        for (xsx, ysy, weight) in diags:
            pairwise.add( (xsx,ysy) )
            pairwise.add( (utils.myGraph.revGene(ysy),utils.myGraph.revGene(xsx)) )
        n_contigs = 0
        n_scaffolds = 0
        n_singletons = 0
        for s in block_names:
            chrom = genome.lstGenes[s]
            # Calling these contigs wouldn't be wrong, but misleading
            if len(chrom) == 1:
                n_singletons += 1
                names[s] = "singleton_%d" % n_singletons
                continue
            # Iterate over each pair of consecutive (stranded) genes
            g1 = chrom[0]
            gs1 = (g1.names[0], g1.strand)
            for g2 in itertools.islice(chrom, 1, None):
                gs2 = (g2.names[0], g2.strand)
                # If one pair is not found, it's a scaffold - that's it !
                if (gs1,gs2) not in pairwise:
                    n_scaffolds += 1
                    names[s] = "scaffold_%d" % n_scaffolds
                    break
                g1 = g2
            else:
                n_contigs += 1
                names[s] = "contig_%d" % n_contigs
    else:
        # Simply use integers
        for (i, s) in enumerate(block_names):
            names[s] = i

    # Print the ancestral genome, using the actual gene names provided by ancGenesFiles
    ancGenes = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc]).lstGenes[None]
    ancGenomeFile = utils.myFile.openFile(arguments["OUT.ancGenomes"] % phylTree.fileName[anc], "w")
    for s in block_names:
        for gene in genome.lstGenes[s]:
            print >> ancGenomeFile, utils.myFile.myTSV.printLine([names[s], gene.beginning, gene.end, gene.strand, " ".join(ancGenes[gene.names[0]].names)])
    ancGenomeFile.close()

start = time.time()
n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)
print >> sys.stderr, "Time elapsed:", time.time() - start
