#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Transform reconstructed ancestral lists of blocks (diags) in formated ancestral genomes (lists of genes)

    usage:
        src/postprocessing/misc.convertContigsToGenome.py example/data/Species.conf A0 \
                -IN.ancDiags=example/results/integrDiags/final/diags.%s.list.bz2 \
                -OUT.ancGenomes=example/results/ancGenomes/final/ancGenome.%s.list.bz2 \
                -ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2
"""

import utils.myGenomes
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs(
    [("phylTree.conf", file), ("target", str)],
    [("IN.ancDiags", str, ""), ("ancGenesFiles", str, ""), ("OUT.ancGenomes", str, "ancGenomes/ancGenome.%s.list.bz2")],
    __doc__
)

# Load species tree - target ancestral genome and the extant species used to assemble blocs
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
targets = phylTree.getTargetsAnc(arguments["target"])

for anc in targets:
    ancGenes = utils.myGenomes.Genome(arguments["ancGenesFiles"] % phylTree.fileName[anc])
    genome = utils.myGenomes.Genome(arguments["IN.ancDiags"] % phylTree.fileName[anc], ancGenes=ancGenes)
    ancGenomeFile = utils.myFile.openFile(arguments["OUT.ancGenomes"] % phylTree.fileName[anc], "w")
    genome.printEnsembl(ancGenomeFile)
    ancGenomeFile.close()
