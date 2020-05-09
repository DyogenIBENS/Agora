#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2016 IBENS/Dyogen : Matthieu MUFFATO, Alexandra Louis, Nga Thi tuy Nguyen, Joseph Lucas, Hugues ROEST CROLLIUS
# mail : hrc@ens.fr or alouis@biologie.ens.fr
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
	Just a script to copy the ancestral blocks of robusts genes (skeleton) in the good directory to continue the procedure
"""

import multiprocessing
import time
import sys

from joblib import Parallel, delayed

import utils.myFile
import utils.myPhylTree
import utils.myTools

# Arguments
arguments = utils.myTools.checkArgs([("phylTree.conf", file), ("target", str)],
                                    [("IN.ancDiags", str, ""), ("OUT.ancDiags", str, ""), ("nbThreads", int, 0)],
                                    __doc__
                                    )

start = time.time()
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
targets = phylTree.getTargetsAnc(arguments["target"])


def do(anc):
    fi = utils.myFile.openFile(arguments["IN.ancDiags"] % phylTree.fileName[anc], "r")
    fo = utils.myFile.openFile(arguments["OUT.ancDiags"] % phylTree.fileName[anc], "w")
    for l in fi:
        print >> fo, l,
    fo.close()
    fi.close()


n_cpu = arguments["nbThreads"] or multiprocessing.cpu_count()
Parallel(n_jobs=n_cpu)(delayed(do)(anc) for anc in targets)
print >> sys.stderr, "Time elapsed:", time.time() - start
