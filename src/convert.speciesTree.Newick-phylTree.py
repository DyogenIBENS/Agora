#! /usr/bin/env python
#  -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """ Convert a species phylogenetic tree from/to newick or phyltree text format

     usage:
        convert.SpeciesTree.Newick-phylTree.py Species.conf -fromNewick > Species.nwk
        convert.SpeciesTree.Newick-phylTree.py Species.nwk +fromNewick > Species.conf """

import utils.myPhylTree as myPhylTree
import utils.myTools as myTools


arguments = myTools.checkArgs([("phylTree.conf",file)], [("fromNewick",bool,True)], __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

if arguments["fromNewick"]:
    phylTree.printPhylTree()
else:
    phylTree.printNewick()
