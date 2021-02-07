#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
        Convert gene trees from NHX to my protTree forest
"""

import sys

import utils.myProteinTree as myProteinTree
import utils.myTools as myTools

arguments = myTools.checkArgs( [("tree",file)], [], __doc__)

for tree in myProteinTree.loadTree(arguments["tree"]):
    myProteinTree.printTree(sys.stdout, tree.data, tree.info, tree.root)
