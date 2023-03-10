#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# AGORA v3.1
# python 3.5
# Copyright © 2006-2022 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021-2022 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import sys

import utils.myAgoraWorkflow

__doc__ = """
    Run the basic AGORA workflow (two single-integration reconstructions)

    Usage:
          src/agora-basic.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2
          src/agora-basic.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -target=A0
          src/agora-basic.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 example/data/genes/genes.%s.list.bz2 -workingDir=example/results -nbThreads=1
"""

(workflow, arguments) = utils.myAgoraWorkflow.AgoraWorkflow.initFromCommandLine(__doc__, [])

workflow.addAncGenesGenerationAnalysis()
workflow.addPairwiseAnalysis(workflow.allAncGenesName)
workflow.addIntegrationAnalysis("denovo", ['+searchLoops'], workflow.allAncGenesName)
workflow.addIntegrationAnalysis("scaffolds", [], None)
workflow.publishGenome(outputName="basic-workflow")

# Launching tasks in multiple threads
#####################################
failed = workflow.tasklist.runAll(arguments["nbThreads"], arguments["sequential"], arguments["forceRerun"])
sys.exit(failed)
