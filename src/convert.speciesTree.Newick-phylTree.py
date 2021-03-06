#! /usr/bin/env python
#  -*- coding: utf-8 -*-
# AGORA v2.1
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """ Convert a species phylogenetic tree from/to newick or phyltree text format

     usage:
        convert.SpeciesTree.Newick-phylTree.py Species.conf -fromNewick > Species.nwk
        convert.SpeciesTree.Newick-phylTree.py Species.nwk +fromNewick > Species.conf """

import utils.myPhylTree as myPhylTree
import utils.myTools as myTools


arguments = myTools.checkArgs([("speciesTree",file)], [("fromNewick",bool,True)], __doc__)

phylTree = myPhylTree.PhylogeneticTree(arguments["speciesTree"])

if arguments["fromNewick"]:
    phylTree.printPhylTree()
else:
    phylTree.printNewick()
