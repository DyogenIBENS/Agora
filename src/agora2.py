#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v1.4
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import multiprocessing
import os
import sys

import utils.myAgoraWorkflow
import utils.myPhylTree
import utils.myTools

__doc__ = """
    Run all the AGORA programs with the v2 preset (with robust families)

    Usage:
          src/agora2.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2
          src/agora2.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -target=A0 -minSize=0.9 -maxSize=1.1
          src/agora2.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -workingDir=example/results -nbThreads=1
"""

arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("geneTrees", file), ("genes", str)],
    [("minSize", float, 1.0), ("maxSize", float, 1.0), ("target", str, ""), ("extantSpeciesFilter", str, ""),
     ("workingDir", str, "."), ("nbThreads", int, multiprocessing.cpu_count())],
    __doc__)

# Path configuration
files = {}
for f in utils.myAgoraWorkflow.AgoraWorkflow.inputParams:
    files[f] = arguments[f]
outputDir = arguments["workingDir"]
for (f, s) in utils.myAgoraWorkflow.AgoraWorkflow.defaultPaths.iteritems():
    files[f] = os.path.normpath(os.path.join(outputDir, s))
scriptDir = os.path.dirname(os.path.abspath(__file__))

dirnameTemplate = "size-%s-%s"
ancGenesDirName = dirnameTemplate % (arguments['minSize'], arguments['maxSize'])

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])
# Check that the syntax is correct
if arguments["target"]:
    phylTree.getTargetsAnc(arguments["target"])
if arguments["extantSpeciesFilter"]:
    phylTree.getTargetsSpec(arguments["extantSpeciesFilter"])

workflow = utils.myAgoraWorkflow.AgoraWorkflow(arguments["target"] or phylTree.root, arguments["extantSpeciesFilter"], scriptDir, files)
workflow.addAncGenesGenerationAnalysis()
workflow.addAncGenesFilterAnalysis(ancGenesDirName, "size", [str(arguments['minSize']), str(arguments['maxSize'])], dirnameTemplate)
workflow.addPairwiseAnalysis(workflow.allAncGenesDirName)
workflow.addPairwiseAnalysis(ancGenesDirName)
workflow.addIntegrationAnalysis("denovo", [], ancGenesDirName)
workflow.addIntegrationAnalysis("refine", [], workflow.allAncGenesDirName)
workflow.addIntegrationAnalysis("extend", ["+onlySingletons"], workflow.allAncGenesDirName)
workflow.addIntegrationAnalysis("halfinsert", [], workflow.allAncGenesDirName)
workflow.addIntegrationAnalysis("groups", [], None)
workflow.addIntegrationAnalysis("publish", [], None, taskName="//")

# Launching tasks in multiple threads
#####################################
failed = workflow.tasklist.runAll(arguments["nbThreads"])
sys.exit(failed)
