#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v3.1
# python 3.5
# Copyright © 2006-2021 IBENS/Dyogen, 2020-2021 EMBL-European Bioinformatics Institute, 2021 Genome Research Ltd : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    When the second argument is a forest of gene trees, read it and extract the gene content of every
    ancestral genome in a separate file.
    When it is a file pattern, assume they are orthology groups (and thus already represent ancestral
    genes), and copy the files over to the location defined by OUT.ancGenesFiles.

    Usage:
        src/ALL.extractGeneFamilies.py example/data/Species.nwk example/data/GeneTreeForest.nhx.bz2 \
                -OUT.ancGenesFiles=example/results/ancGenes/all/ancGenes.%s.list.bz2
"""

import collections
import os
import sys

import utils.myFile
import utils.myPhylTree
import utils.myProteinTree
import utils.myTools

sys.setrecursionlimit(10000)

# Arguments
###########

arguments = utils.myTools.checkArgs(
    [("speciesTree", utils.myTools.FileArgChecker), ("geneTrees|orthologyGroups", utils.myTools.FileOrPatternArgChecker)],
    [("OUT.ancGenesFiles", str, ""), ("reuseNames", bool, False)],
    __doc__
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

if "%s" in arguments["geneTrees|orthologyGroups"]:
    # ancGenes are already present, just copy them over, but add a family name for unique reference
    for anc in phylTree.listAncestr:
        inputPath = arguments["geneTrees|orthologyGroups"] % phylTree.fileName[anc]
        if os.path.exists(inputPath):
            if anc in phylTree.listSpecies:
                print("Copying families of genome %s ..." % anc, end=' ', file=sys.stderr)
            else:
                print("Copying families of ancestral genome %s ..." % anc, end=' ', file=sys.stderr)
            code = 'FAM' + anc[:4].upper()
            outputPath = arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc]
            fi = utils.myFile.openFile(inputPath, "r")
            fo = utils.myFile.openFile(outputPath, "w")
            n = 0
            for l in fi:
                n += 1
                fo.write(code + ('%06d ' % n) + l)
            fo.close()
            fi.close()
            print(n, "OK", file=sys.stderr)
        else:
            print("No file for '%s' in '%s'" % (anc, arguments["geneTrees|orthologyGroups"]))
    sys.exit(0)


dupCount = collections.defaultdict(int)


def futureName(name, dup):
    if dup >= 2:
        dupCount[name] += 1
        return name + utils.myProteinTree.getDupSuffix(dupCount[name], False)
    else:
        return name


#########################
# Find roots in families #
#########################
def getRoots(node, previousAnc, lastWrittenAnc):
    newAnc = tree.info[node]['taxon_name']
    (_, newLastWritten, isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc, newAnc,
                                                                         tree.info[node]['Duplication'] >= 2)

    if isroot:
        return [node]

    # descendant genes
    subRoots = []
    for (g, _) in tree.data.get(node, []):
        subRoots.extend(getRoots(g, newAnc, newLastWritten))
    return subRoots


count = collections.defaultdict(int)


#################################
# Backup all the gene families #
################################
def extractGeneFamilies(node, baseName, previousAnc, lastWrittenAnc):
    newAnc = tree.info[node]['taxon_name']
    (toWrite, newLastWritten, isroot) = utils.myProteinTree.getIntermediateAnc(phylTree, previousAnc, lastWrittenAnc,
                                                                               newAnc,
                                                                               tree.info[node]['Duplication'] >= 2)

    if isroot and (previousAnc is not None):
        if not arguments["reuseNames"]:
            baseName = baseName.split(".")[0]
        count[baseName] += 1
        currName = baseName + utils.myProteinTree.getDupSuffix(count[baseName], True)
    else:
        currName = baseName
    tree.info[node]['family_name'] = currName

    # descendant genes
    if node in tree.data:
        allGenes = []
        for (g, _) in tree.data[node]:
            allGenes.extend(
                extractGeneFamilies(g, futureName(currName, tree.info[node]['Duplication']), newAnc, newLastWritten))
    else:
        allGenes = [tree.info[node]["gene_name"]]

    for a in toWrite:
        geneFamilies[a].append([currName] + allGenes)

    return allGenes


geneFamilies = collections.defaultdict(list)
for tree in utils.myProteinTree.loadTree(arguments["geneTrees|orthologyGroups"]):
    extractGeneFamilies(tree.root, tree.info[tree.root]["tree_name"], None, None)
    if tree.info[tree.root]["format"] == "NHX":
        tree.printNewick(sys.stdout, withDist=True, withTags=True, withAncSpeciesNames=True, withAncGenesNames=True)
    else:
        tree.printTree(sys.stdout)

for anc in sorted(geneFamilies):
    if anc in phylTree.listSpecies:
        print("Writing families of genome %s ..." % anc, end=' ', file=sys.stderr)
    else:
        print("Writing families of ancestral genome %s ..." % anc, end=' ', file=sys.stderr)
    lst = geneFamilies[anc]
    f = utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc], "w")
    for gg in lst:
        print(" ".join(gg), file=f)
    f.close()
    print(len(lst), "OK", file=sys.stderr)
