#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.0
# python 2.7
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import multiprocessing
import os
import sys

import utils.myAgoraWorkflow
import utils.myPhylTree
import utils.myTools
from utils.myTools import file

__doc__ = """
    Run the basic AGORA workflow (two single-integration reconstructions)

    Usage:
          src/agora-basic.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2
          src/agora-basic.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -target=A0
          src/agora-basic.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -workingDir=example/results -nbThreads=1
"""

arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("geneTrees", file), ("genes", str)],
    [("target", str, ""), ("extantSpeciesFilter", str, ""),
     ("workingDir", str, "."), ("nbThreads", int, multiprocessing.cpu_count()), ("forceRerun", bool, False), ("sequential", bool, True)],
    __doc__)

# Path configuration
files = {}
for f in utils.myAgoraWorkflow.AgoraWorkflow.inputParams:
    files[f] = arguments[f]
outputDir = arguments["workingDir"]
for (f, s) in utils.myAgoraWorkflow.AgoraWorkflow.defaultPaths.items():
    files[f] = os.path.normpath(os.path.join(outputDir, s))
scriptDir = os.path.dirname(os.path.abspath(__file__))

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
# Check that the syntax is correct
if arguments["target"]:
    phylTree.getTargetsAnc(arguments["target"])
if arguments["extantSpeciesFilter"]:
    phylTree.getTargetsSpec(arguments["extantSpeciesFilter"])

workflow = utils.myAgoraWorkflow.AgoraWorkflow(arguments["target"] or phylTree.root, arguments["extantSpeciesFilter"], scriptDir, files)
workflow.addAncGenesGenerationAnalysis()
workflow.addPairwiseAnalysis(workflow.allAncGenesName)
workflow.addIntegrationAnalysis("denovo", ['+searchLoops'], workflow.allAncGenesName)
workflow.addIntegrationAnalysis("scaffolds", [], None)
workflow.publishGenome(outputName="basic-workflow")

# Launching tasks in multiple threads
#####################################
failed = workflow.tasklist.runAll(arguments["nbThreads"], arguments["sequential"], arguments["forceRerun"])
sys.exit(failed)
