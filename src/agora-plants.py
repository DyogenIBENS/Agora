#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.1
# python 3.5
# Copyright © 2006-2022 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021-2022 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import sys

import utils.myAgoraWorkflow

__doc__ = """
    Run the Plants workflow of AGORA (two multi-integration passes)

    Usage:
          src/agora-plants.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2
          src/agora-plants.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -target=A0 -minSize=0.9 -maxSize=1.1
          src/agora-plants.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -workingDir=example/results -nbThreads=1
"""

options = [("minSize", float, 0.9), ("maxSize", float, 1.1)]
(workflow, arguments) = utils.myAgoraWorkflow.AgoraWorkflow.initFromCommandLine(__doc__, options)

workflow.addAncGenesGenerationAnalysis()
workflow.reconstructionPassWithAncGenesFiltering("size", [arguments['minSize'], arguments['maxSize']])

workflow.useBlocksAsAncGenes()

workflow.addPairwiseAnalysis(workflow.allAncGenesName, params=["-anchorSize=3"])
workflow.addIntegrationAnalysis("denovo", [], workflow.allAncGenesName)
workflow.convertToRealAncGenes()
workflow.markForSelection()

filtBlocksMethods = [("propLength", "50"), ("propLength", "70"), ("fixedLength", "20"), ("fixedLength", "50")]
for filtParams in filtBlocksMethods:
    workflow.reconstructionPassWithAncGenesFiltering(filtParams[0], list(filtParams[1:]))
    workflow.convertToRealAncGenes()
    workflow.markForSelection()

workflow.revertToRealAncGenes()
workflow.addSelectionAnalysis(taskName="best-pass2")
workflow.publishGenome(outputName="plants-workflow")

# Launching tasks in multiple threads
#####################################
failed = workflow.tasklist.runAll(arguments["nbThreads"], arguments["sequential"], arguments["forceRerun"])
sys.exit(failed)
