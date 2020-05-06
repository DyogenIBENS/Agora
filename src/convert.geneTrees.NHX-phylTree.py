#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
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