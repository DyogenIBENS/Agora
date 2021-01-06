#! /usr/bin/env python
# -*- coding: utf-8 -*-
# AGORA v2.0
# python 2.7
# Copyright © 2021 IBENS/Dyogen and EMBL-European Bioinformatics Institute : Matthieu MUFFATO, Alexandra LOUIS, Thi Thuy Nga NGUYEN, Hugues ROEST CROLLIUS
# mail : agora@bio.ens.psl.eu
# This is free software; you may copy, modify and/or distribute this work under the terms of the GNU General Public License, version 3 or later and the CeCiLL v2 license in France

__doc__ = """
    Read the forest of gene trees and extract each ancestral genome gene content in separate file.
    Extract ancestral gene content from a forest of gene trees.
    One file per ancestor, and one file per extant species.

        Usage:
            ALL.extractGeneFamilies.py PhylTree.conf GeneTrees.bz
  """

import collections
import sys

import utils.myFile
import utils.myPhylTree
import utils.myProteinTree
import utils.myTools

sys.setrecursionlimit(10000)

# Arguments
###########

arguments = utils.myTools.checkArgs(
    [("speciesTree", file), ("geneTrees", file)],
    [("OUT.ancGenesFiles", str, ""), ("reuseNames", bool, False)],
    __doc__
)

phylTree = utils.myPhylTree.PhylogeneticTree(arguments["speciesTree"])

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
for tree in utils.myProteinTree.loadTree(arguments["geneTrees"]):
    extractGeneFamilies(tree.root, tree.info[tree.root]["tree_name"], None, None)
    if tree.info[tree.root]["format"] == "NHX":
        tree.printNewick(sys.stdout, withDist=True, withTags=True, withAncSpeciesNames=True, withAncGenesNames=True)
    else:
        tree.printTree(sys.stdout)

for (anc, lst) in geneFamilies.iteritems():
    if anc in phylTree.listSpecies:
        print >> sys.stderr, "Writing families of genome %s ..." % anc,
    else:
        print >> sys.stderr, "Writing families of ancestral genome %s ..." % anc,
    f = utils.myFile.openFile(arguments["OUT.ancGenesFiles"] % phylTree.fileName[anc], "w")
    for gg in lst:
        print >> f, " ".join(gg)
    f.close()
    print >> sys.stderr, len(lst), "OK"
