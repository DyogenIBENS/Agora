#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.0
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import multiprocessing
import os
import sys

import utils.myAgoraWorkflow
import utils.myPhylTree
import utils.myTools

__doc__ = """
    Run the Vertebrates workflow of AGORA (i.e. with constrained families)

    Usage:
          src/agora-plants.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2
          src/agora-plants.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -target=A0 -minSize=0.9 -maxSize=1.1
          src/agora-plants.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -workingDir=example/results -nbThreads=1
"""

arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("geneTrees", file), ("genes", str)],
    [("minSize", float, 1.0), ("maxSize", float, 1.0), ("target", str, ""), ("extantSpeciesFilter", str, ""),
     ("workingDir", str, "."), ("nbThreads", int, multiprocessing.cpu_count()), ("sequential", bool, True)],
    __doc__)

# Path configuration
files = {}
for f in utils.myAgoraWorkflow.AgoraWorkflow.inputParams:
    files[f] = arguments[f]
outputDir = arguments["workingDir"]
for (f, s) in utils.myAgoraWorkflow.AgoraWorkflow.defaultPaths.iteritems():
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

workflow.reconstructionPassWithAncGenesFiltering("size", [arguments['minSize'], arguments['maxSize']])
# workflow.addIntegrationAnalysis("scaffolds", [], None)
workflow.useBlocksAsAncGenes()

# TODO denovo all and compare to "scaffolds"
workflow.addPairwiseAnalysis(workflow.allAncGenesName, params=["-anchorSize=3"])
workflow.addIntegrationAnalysis("denovo", [], workflow.allAncGenesName)
workflow.convertToRealAncGenes()
workflow.markForSelection()

filtBlocksMethods = [("propLength", "50"), ("propLength", "70"), ("fixedLength", "20"), ("fixedLength", "50")]

for filtParams in filtBlocksMethods:
    workflow.reconstructionPassWithAncGenesFiltering(filtParams[0], list(filtParams[1:]))
    workflow.convertToRealAncGenes()
    workflow.markForSelection()
# TODO name the output
workflow.addSelectionAnalysis()
# workflow.publishGenome(outputName="plants-workflow")

# workflow.tasklist.printGraphviz(sys.stdout)
# sys.exit(0)
# Launching tasks in multiple threads
#####################################
failed = workflow.tasklist.runAll(arguments["nbThreads"], arguments["sequential"])
sys.exit(failed)
