#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Transform a reconstructed ancestral list of blocks (diags) in a formated ancestral genome (list of genes)

    usage:
        ./misc.convertContigsToGenome.py diags.Amiota.list.bz2 ancGenes.Amniota.list.bz2 > genome.Amniota.list
"""

import sys

import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs([("contigsFile", file), ("ancGenesFile", file)], [], __doc__)

ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"])

genome = utils.myGenomes.Genome(arguments["contigsFile"], ancGenes=ancGenes)

genome.printEnsembl(sys.stdout)
