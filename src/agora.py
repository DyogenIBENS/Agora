#! /usr/bin/env python
# -*- coding: utf-8 -*-
# Agora v1.00
# python 2.7
# Copyright © 2020 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

import multiprocessing
import os
import re
import subprocess
import sys

import utils.myAgoraWorkflow
import utils.myFile
import utils.myPhylTree
import utils.myTools

__doc__ = """
    Run all the AGORA programs thanks to configuration file

    Usage:
          src/agora.py conf/agora.ini
          src/agora.py conf/agora.ini -workingDir=example/results -nbThreads=4
"""

arguments = utils.myTools.checkArgs(
    [("agora.conf", file)],
    [("workingDir", str, "."), ("nbThreads", int, multiprocessing.cpu_count()),
     ],
    __doc__)

# loading configuration file
################################
bysections = {}
f = utils.myFile.openFile(arguments["agora.conf"], "r")
for l in f:

    l = l.partition("#")[0].strip()
    if l.startswith(">"):
        curr = []
        bysections[l[1:].strip().lower()] = curr
    elif len(l) > 1:
        curr.append(l)
f.close()


# Cut a string for delim
##########################
def partition(s, delim):
    x = s.partition(delim)
    return (x[0].strip(), x[2].strip())


# Files and directory section
##############################
files = {}
for x in bysections["files"]:
    x = partition(x, "=")
    files[x[0].lower()] = x[1]

# All input paths are relative to the directory of the configuration file
inputDir = os.path.dirname(arguments["agora.conf"])
outputDir = arguments["workingDir"]
inputParams = ["speciestree", "genes", "genetrees"]
for f in files:
    files[f] = os.path.normpath(os.path.join(inputDir if f in inputParams else outputDir, files[f]))
scriptDir = os.path.dirname(os.path.abspath(__file__))

phylTree = utils.myPhylTree.PhylogeneticTree(files["speciestree"])


workflow = utils.myAgoraWorkflow.AgoraWorkflow(phylTree.root, scriptDir, files)

# Ancestral genes lists Section
################################


# all ancGenes task - gather nickname and only launch if explicitly requested to (backwards compatibility)
allname = "all"
launchall = False
patternall = re.compile(r'=.*\ball\b')
for x in bysections["ancgenes"]:
    if patternall.search(x):
        x = partition(x, "=")
        launchall = "*" not in x[0]
        allname = x[0].replace("*", "").strip()

workflow.addAncGenesGenerationAnalysis(launchall)

ancGenes = {allname: "all", "0": "all"}

# Parse the section
for x in bysections["ancgenes"]:
    if patternall.search(x):
        continue
    x = partition(x, "=")
    (params, root) = partition(x[1], "!")
    if root == "":
        root = phylTree.root
    # index
    t = x[0].replace("*", "").split(",")
    sizes = params.split(",")
    assert len(t) == len(sizes)
    minSizes = []
    maxSizes = []
    dirnameTemplate = "size-%s-%s"
    ancGenesDirNames = []
    for i in range(len(sizes)):
        size = sizes[i].split()
        minSizes.append(size[0])
        maxSizes.append(size[1])
        dirname = dirnameTemplate % tuple(size)
        ancGenes[t[i]] = dirname
        ancGenesDirNames.append(dirname)
    minSizesStr = ",".join(minSizes)
    maxSizesStr = ",".join(maxSizes)

    taskname = dirnameTemplate % (minSizesStr, maxSizesStr)
    workflow.addAncGenesFilterAnalysis(taskname, "size", dirnameTemplate, root, [minSizesStr, maxSizesStr], "*" not in x[0])

    if len(ancGenesDirNames) > 1:
        for dirname in ancGenesDirNames:
            workflow.addDummy(("ancgenes", dirname), [("ancgenes", taskname)])

# Pairwise comparison section
#############################


for x in bysections.get("pairwise", []):
    x = partition(x, "=")
    dirname = ancGenes[x[0].replace("*", "").strip()]
    (params, root) = partition(x[1], "!")
    params = params.split()
    if root == "":
        root = phylTree.root

    # Pairwise comparison tasks
    workflow.addPairwiseAnalysis(dirname, params[0], root, params[1:], "*" not in x[0])

# Integration section
#####################
for x in bysections.get("integration", []):
    tolaunch = ("*" not in x)

    (params, root) = partition(x, "!")
    (params, input) = partition(params, "<")
    (params, output) = partition(params, ">")
    params = params.replace("*", "").split()

    currMethod = None
    # method's name
    if params[0].startswith("["):
        currMethod = params.pop(0)[1:-1]

    # used compared pairs
    dirname = None
    if params[-1].startswith("("):
        dirname = ancGenes[params.pop()[1:-1]]

    workflow.addIntegrationAnalysis(currMethod, params[0], input, output, dirname, root, params[1:], tolaunch)

# Launching tasks in multiple threads
#####################################
failed = workflow.tasklist.runAll(arguments["nbThreads"])
sys.exit(failed)
