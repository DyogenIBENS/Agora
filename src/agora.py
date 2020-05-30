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


tasklist = utils.myAgoraWorkflow.TaskList()

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

tasklist.addTask(
    ("ancgenes", "all"),
    [],
    (
        [os.path.join(scriptDir, "ALL.extractGeneFamilies.py"), files["speciestree"], files["genetrees"],
            "-OUT.ancGenesFiles=" + files["ancgenesdata"] % {"filt": "all", "name": "%s"}],
        files["genetreeswithancnames"],
        files["ancgeneslog"] % {"filt": "ancGenes"},
        launchall,
    )
)

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
    ancGenesDirNames = []
    for i in range(len(sizes)):
        size = sizes[i].split()
        minSizes.append(size[0])
        maxSizes.append(size[1])
        dirname = "size-" + str(size[0]) + "-" + str(size[1])
        ancGenes[t[i]] = dirname
        ancGenesDirNames.append(dirname)
    minSizesStr = ",".join(minSizes)
    maxSizesStr = ",".join(maxSizes)

    taskname = "size-" + minsize + "-" + maxsize
    tasklist.addTask(
        ("ancgenes", taskname),
        [("ancgenes", "all")],
        (
            [os.path.join(scriptDir, "ALL.filterGeneFamilies-%s.py" % "size"), files["speciestree"], root,
             files["ancgenesdata"] % {"filt": "all", "name": "%s"},
             files["ancgenesdata"] % {"filt": "size-%s-%s", "name": "%s"}] + [minSizesStr, maxSizesStr],
            os.devnull,
            files["ancgeneslog"] % {"filt": taskname},
            "*" not in x[0]
        )
    )

    if len(ancGenesDirNames) > 1:
        for dirname in ancGenesDirNames:
            tasklist.addTask( ("ancgenes", dirname), [("ancgenes", taskname)], (None, None, None, False) )

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
    tasklist.addTask(
        ("pairwise", dirname),
        [("ancgenes",dirname)],
        (
            [os.path.join(scriptDir, "buildSynteny.pairwise-%s.py" % params[0]), files["speciestree"], root,
             "-ancGenesFiles=" + files["ancgenesdata"] % {"filt": dirname, "name": "%s"},
             "-genesFiles=" + files["genes"] % {"name": "%s"},
             "-OUT.pairwise=" + files["pairwiseoutput"] % {"filt": dirname, "name": "%s"}] + params[1:],
            os.devnull,
            files["pairwiselog"] % {"filt": dirname},
            "*" not in x[0]
        )
    )

# Integration section
#####################
interm = {}
refMethod = {}
for x in bysections.get("integration", []):
    tolaunch = ("*" not in x)

    (params, root) = partition(x, "!")
    (params, input) = partition(params, "<")
    (params, output) = partition(params, ">")
    params = params.replace("*", "").split()

    # method's name
    currMethod = params.pop(0)[1:-1] if params[0].startswith("[") else params[0]

    # used compared pairs
    dirname = None
    if params[-1].startswith("("):
        dirname = ancGenes[params.pop()[1:-1]]
        currMethod = currMethod[:-1] if currMethod.endswith("/") else (currMethod + "-" + dirname)
    elif currMethod.endswith("/"):
        currMethod = currMethod[:-1]

    if params[0] == "denovo":
        newMethod = currMethod
        if root == "":
            root = phylTree.root
        refMethod[newMethod] = (newMethod, root)
    else:
        # input data
        prevMethod = interm[input] if input != "" else newMethod
        newMethod = currMethod[1:] if currMethod.startswith("/") else (prevMethod + "." + currMethod)
        if root == "":
            root = refMethod[prevMethod][1]
        refMethod[newMethod] = (newMethod, root) if params[0] == "refine" else refMethod[prevMethod]

    if output != "":
        interm[output] = newMethod

    # task parameters
    args = [os.path.join(scriptDir, "buildSynteny.integr-%s.py" % params[0]), files["speciestree"], root] + params[1:]

    if params[0] == "publish":
        # "publish" is not an integration method
        args[0] = os.path.join(scriptDir, "convert.ancGenomes.diags-genes.py")
        args.append("-OUT.ancGenomes=" + files["ancgenomesoutput"] % {"method": newMethod, "name": "%s"})
        logfile = files["ancgenomeslog"]
    else:
        args.append("-OUT.ancDiags=" + files["integrblocks"] % {"method": newMethod, "name": "%s"})
        logfile = files["integrlog"]

    dep = []
    if dirname is not None:
        dep.append(("pairwise", dirname))
        args.append(files["pairwiseoutput"] % {"filt": dirname, "name": "%s"})

    if params[0] in ["denovo", "groups", "publish"]:
        args.append("-ancGenesFiles=" + files["ancgenesdata"] % {"filt": "all", "name": "%s"})

    # No input data to consider for the denovo method
    if params[0] != "denovo":
        dep.append(("integr", prevMethod))
        args.append("-IN.ancDiags=" + files["integrblocks"] % {"method": prevMethod, "name": "%s"})

    if params[0] == "halfinsert":
        # The script needs singleton reference for "halfinsert"
        dep.append(("integr", refMethod[newMethod][0]))
        args.append("-REF.ancDiags=" + files["integrblocks"] % {"method": refMethod[newMethod][0], "name": "%s"})

    if params[0] == "groups":
        args.append("-genesFiles=" + files["genes"] % {"name": "%s"})

    if params[0] not in ["copy", "publish"]:
        args.append("-LOG.ancGraph=" + files["integroutput"] % {"method": newMethod, "name": "%s"})

    # Currently, all those methods are multithreaded
    multithreaded = True

    # Integration task
    tasklist.addTask(
        ("integr", newMethod),
        dep,
        (
            args,
            os.devnull,
            logfile % {"method": newMethod},
            tolaunch
        ),
        multithreaded
    )

    # The publish method doesn't generate integrDiags and can't be used as an input method
    if params[0] == "publish":
        newMethod = prevMethod

# Launching tasks in multiple threads
#####################################
failed = tasklist.runAll(arguments["nbThreads"])
sys.exit(failed)
